#include "pthreads.hpp"
#include "cstdlib"
#include <iostream>
#include <memory>

inline void numClusters(size_t &nclus, const std::vector<size_t> &c)
{
    size_t n = c.size();
    std::vector<size_t> discreetClus(n, 0); // vector where the ith element is a if cluster i has a nodes

    for (size_t i = 0; i < n; i++)
    {
        discreetClus[c[i] - 1] = 1; // we assume that there is no ith row and column that are both zero so we know that all ids included in c exist in A
                                    // we can atomically add 1, instead, to the cluster of the ith row to know how many nodes are in each cluster
    }

    nclus = 0;
    for (size_t i = 0; i < n; i++)
    {
        if (discreetClus[i] == 0)
            continue;
        nclus += 1;
    }
}

inline void *fnClearAux(void *args)
{
    forThread *forArgs = (forThread *)args;
    for (size_t i = forArgs->start; i < forArgs->end; i++)
        forArgs->array[i] = 0; // reset auxiliary vector

    return nullptr;
}

void *fnSumAux(void *args)
{
    sumThread *sumArgs = (sumThread *)args;
    size_t n = sumArgs->row.size();

    size_t end = sumArgs->end;
    if (end > n)
        end = n;

    size_t x, endAux;
    for (size_t i = sumArgs->start; i < end; i++)
    {
        if (sumArgs->id != sumArgs->c[i]) // c[i] : cluster of row i of colCompressed / row
            continue;

        if (i == (n - 1)) // last row contains the last range of non-zero elements
            endAux = sumArgs->col.size();
        else
            endAux = sumArgs->row[i + 1];

        x = sumArgs->row[i];
        if (x == endAux) // this row has no non-zero elements
            continue;

        if (!sumArgs->clusterHasElements) // declare that this row has non-zero elements only once
            sumArgs->clusterHasElements = 1;

        for (size_t j = x; j < endAux; j++)
        {
            sumArgs->auxValueVector[sumArgs->c[sumArgs->col[j]] - 1].fetch_add(sumArgs->val[j]); // compress cols by summing the values of each cluster to the column the cluster id points to
        }
    }

    return nullptr;
}

void *fnAssignM(void *args)
{

    size_t localCount;
    assignThread *assignArgs = (assignThread *)args;
    for (size_t i = assignArgs->start; i < assignArgs->end; i++)
    {
        if (assignArgs->auxValueVector[i] == 0)
            continue;

        localCount = assignArgs->allCount->fetch_add(1); // returns the previous value

        assignArgs->valM[localCount] = assignArgs->auxValueVector[i];
        assignArgs->colM[localCount] = i;
    }

    return nullptr;
}

void pthreads(std::vector<size_t> &rowM, std::vector<size_t> &colM, std::vector<uint32_t> &valM,
              const std::vector<size_t> &row, const std::vector<size_t> &col, const std::vector<uint32_t> &val, const std::vector<size_t> &c)
{
    if (row.size() != c.size())
    {
        printf("Error: sizes of row and c are incompatible\n");
        exit(1);
    }
    else if (col.size() != val.size())
    {
        printf("Error: sizes of col and val are incompatible\n");
        exit(1);
    }
    else if (row.size() == col.size())
    {
        printf("Error: CSR requires more space than dense matrix representation \n Use dense matrix implementation instead...\n");
        exit(1);
    }
    else if ((rowM.size() != row.size()) || (colM.size() != col.size()) || (valM.size() != val.size()))
    // if the input graph is its the minor, the dimensions of the compressed vectors should be equal to the dimensions of the input vectors
    {
        printf("Error: at least one of the compressed vectors doesn't have enough allocated space \n");
        exit(1);
    }

    size_t n = c.size();
    size_t nz = val.size();
    size_t nclus;

    numClusters(nclus, c);

    std::atomic<size_t> allCount(0);
    bool clusterHasElements = 0;
    std::vector<std::atomic<uint32_t>> auxValueVector(nclus); // auxiliary vector that will contain all the non-zero values of each cluster (element of rowM)
    rowM.resize(nclus);                                       // resize vector to the number of clusters

    size_t cacheLines = n / 16; // how many cache lines the array fills:
                                // a chunk of x elements of an int array (4*x bytes) will be equal to a cache line (64 bytes) to avoid false sharing
    size_t numThreads = 4;      // consider edge case where we need + 1 thread

    if (!cacheLines)
    {
        numThreads = 1;
        cacheLines = 1;
    }
    else if (cacheLines <= numThreads) // the array of the smallest type (int < size_t) fits in <=4 cache lines
    {
        numThreads = cacheLines;
    }

    size_t chunk = cacheLines * ELEMENTS_PER_CACHE_LINE / numThreads;
    size_t lastThreadChunk = chunk;
    printf("chunk: %ld\n", chunk);

    if (n > (cacheLines * ELEMENTS_PER_CACHE_LINE)) // n does not fill in an integer number of cachelines
    {
        lastThreadChunk = chunk + (n - cacheLines * ELEMENTS_PER_CACHE_LINE); // the last thread gets its chunk + the remaining items
    }
    else if (n <= (cacheLines * ELEMENTS_PER_CACHE_LINE)) // the array fits in a single cache line
        lastThreadChunk = n;

    printf("last chunk: %ld\n", lastThreadChunk);

    size_t cacheLinesClus = nclus / 16;
    size_t numThreadsClus = 4;
    if (!cacheLinesClus) // the auxValueVector array fits in a cache line
    {
        numThreadsClus = 1;
        cacheLinesClus = 1;
    }
    else if (cacheLinesClus <= numThreadsClus)
    {
        numThreadsClus = cacheLinesClus;
    }

    size_t chunkClus = cacheLinesClus * ELEMENTS_PER_CACHE_LINE / numThreadsClus;
    size_t lastThreadChunkClus = chunkClus;

    if (nclus > (cacheLinesClus * ELEMENTS_PER_CACHE_LINE)) // n does not fill in an integer number of cachelines
    {
        lastThreadChunkClus = chunkClus + (nclus - cacheLinesClus * ELEMENTS_PER_CACHE_LINE); // the last thread gets its chunk + the remaining items
    }
    else if (nclus <= (cacheLinesClus * ELEMENTS_PER_CACHE_LINE)) // the array fits in a single cache line
        lastThreadChunkClus = nclus;

    // if (!numThreadsClus)
    // {
    //     numThreadsClus = 1;
    // }
    // else if (numThreadsClus > 4)
    // {
    //     numThreadsClus = 4;
    // }

    // size_t lastThreadClusChunk = chunksClusPerThread;

    // if (chunkClus > (chunksClusPerThread * 4))
    // {
    //     lastThreadClusChunk = chunkClus - (3 * chunksClusPerThread);
    // }

    std::cout << "Running with numThreads, numThreadsClus: " << numThreads << ", " << numThreadsClus << std::endl;

    pthread_t threads[numThreads], threadsClus[numThreadsClus], threadsClus1[numThreadsClus];

    for (size_t id = 1; id < (nclus + 1); id++) // cluster ids start from 1
    {
        rowM[id - 1] = allCount;
        clusterHasElements = 0;

        std::vector<forThread> forArgs;
        forArgs.reserve(numThreadsClus);

        for (size_t i = 0; i < numThreadsClus; i++)
        {
            size_t start = i * chunkClus;
            size_t end = (i == (numThreadsClus - 1)) ? (start + lastThreadChunkClus) : (start + chunkClus);
            forArgs.push_back({start, end, auxValueVector});
            pthread_create(&threadsClus[i], NULL, fnClearAux, (void *)&forArgs[i]);
        }

        for (pthread_t i : threadsClus)
        {
            pthread_join(i, NULL);
        }

        std::vector<sumThread> sumArgs;
        sumArgs.reserve(numThreads);

        for (size_t i = 0; i < numThreads; i++)
        {
            size_t start = i * chunk;
            size_t end = (i == (numThreads - 1)) ? (start + lastThreadChunk) : (start + chunk);
            sumArgs.push_back({id, start, end, row, col, val, c, auxValueVector, clusterHasElements});
            pthread_create(&threads[i], NULL, fnSumAux, (void *)&sumArgs[i]);
        }

        for (pthread_t i : threads)
        {
            pthread_join(i, NULL);
        }

        if (!clusterHasElements)
            continue;

        std::vector<assignThread> assignArgs;
        assignArgs.reserve(numThreadsClus);

        for (size_t i = 0; i < numThreadsClus; i++)
        {
            size_t start = i * chunk;
            size_t end = (i == (numThreadsClus - 1)) ? (start + lastThreadChunkClus) : (start + chunkClus);
            assignArgs.push_back({start, end, auxValueVector, &allCount, colM, valM});
            pthread_create(&threadsClus1[i], NULL, fnAssignM, (void *)&assignArgs[i]);
        }

        for (pthread_t i : threadsClus1)
        {
            pthread_join(i, NULL);
        }
    }

    colM.resize(allCount);
    valM.resize(allCount);
}