#include "pthreads.hpp"
#include "cstdlib"
#include "iostream"
#include "sys/time.h"
#include <chrono>

void *fnNumClusters(void *args)
{
    nclusThread *nclusArgs = (nclusThread *)args;
    size_t n = nclusArgs->c.size();
    std::vector<size_t> discreetClus(n, 0); // vector where the ith element is a if cluster i has a nodes

    for (size_t i = nclusArgs->start; i < nclusArgs->end; i++)
    {
        discreetClus[(nclusArgs->c[i] - 1)] = 1; // we assume that there is no ith row and column that are both zero so we know that all ids included in c exist in A
                                                 // we can atomically add 1, instead, to the cluster of the ith row to know how many nodes are in each cluster
    }

    for (size_t i = nclusArgs->start; i < nclusArgs->end; i++)
    {
        if (discreetClus[i] == 0)
            continue;
        nclusArgs->nclus++;
    }

    return nullptr;
}

void *fnClearAux(void *args)
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
    size_t index = sumArgs->auxValueIndex;

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
                        sumArgs->auxValueVector[index + sumArgs->c[sumArgs->col[j]] - 1] += (sumArgs->val[j]); // compress cols by summing the values of each cluster to the column the cluster id points to
        }
    }

    return nullptr;
}

void *fnReduction(void *args)
{ // numThreadsClus
    reductionThread *redArgs = (reductionThread *)args;
    size_t nclus = redArgs->nclus;
    size_t size = redArgs->inputArray.size();

    for (size_t id = redArgs->id_start; id < redArgs->id_end; id++)
    {
        for (size_t i = (id - 1); i < size; i += nclus)
        {
            redArgs->outputArray[id - 1] += redArgs->inputArray[i];
        }
    }

    return nullptr;

    // else per nclus*numThreads chunk: redArgs->outputArray[id] += redArgs->inputArray[i];
}

void *fnAssignM(void *args)
{

    size_t localCount;
    assignThread *assignArgs = (assignThread *)args;
    for (size_t i = assignArgs->start; i < assignArgs->end; i++)
    {
        if (assignArgs->auxValueVector[i] == 0)
            continue;

        localCount = assignArgs->allCount.fetch_add(1); // returns the previous value

        assignArgs->valM[localCount] = assignArgs->auxValueVector[i];
        assignArgs->colM[localCount] = i;
    }
    // printf("localCount: %ld\n", localCount);

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

    size_t cacheLines = n / ELEMENTS_PER_CACHE_LINE; // how many cache lines the array fills:
                                                     // a chunk of x elements of an int array (4*x bytes) will be equal to a cache line (64 bytes) to avoid false sharing
    size_t numThreads = 4;                           // consider edge case where we need + 1 thread

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

    if (n > (cacheLines * ELEMENTS_PER_CACHE_LINE)) // n does not fill in an integer number of cachelines
    {
        lastThreadChunk = chunk + (n - cacheLines * ELEMENTS_PER_CACHE_LINE); // the last thread gets its chunk + the remaining items
    }
    else // the array fits in a single cache line
        lastThreadChunk = n;

    // if we used different variables for each thread and then calculate their sum,
    // since this vector of variables would fit in a single cache line,
    // the whole vector would be updated each time an element is changed, so it would be as or more costly

    size_t nclus = 0;
    std::vector<size_t> nclusVector(numThreads, 0);
    std::vector<nclusThread> nclusArgs;
    nclusArgs.reserve(numThreads);

    pthread_t nclusThreads[numThreads];

    for (size_t i = 0; i < numThreads; i++)
    {
        size_t start = i * chunk;
        size_t end = (i == (numThreads - 1)) ? (start + lastThreadChunk) : (start + chunk);
        nclusArgs.push_back({start, end, c, nclusVector[i]});
        pthread_create(&nclusThreads[i], NULL, fnNumClusters, (void *)&nclusArgs[i]);
    }

    for (pthread_t i : nclusThreads)
    {
        pthread_join(i, NULL);
    }

    for (size_t i = 0; i < numThreads; i++)
    {
        nclus += nclusVector[i];
    }

    std::atomic<size_t> allCount(0);
    bool clusterHasElements = 0;
    rowM.resize(nclus); // resize vector to the number of clusters

    size_t cacheLinesClus = nclus / ELEMENTS_PER_CACHE_LINE;
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

    // std::cout << "Running with numThreads, numThreadsClus: " << numThreads << ", " << numThreadsClus << std::endl;

    pthread_t threadsSum[numThreads], threadsAssign[numThreadsClus], reductionThreads[numThreadsClus];

    pthread_attr_t attr;
    cpu_set_t cpu_set;

    pthread_attr_init(&attr);

    // Set CPU affinity to use only CPU 0
    CPU_ZERO(&cpu_set);
    CPU_SET(0, &cpu_set);

    if (pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &cpu_set) != 0)
    {
        perror("pthread_attr_setaffinity_np");
        exit(EXIT_FAILURE);
    }

    std::vector<uint32_t> valMVector(nclus, 0);
    std::vector<uint32_t> auxValueVector(nclus * numThreads, 0);
    // std::vector<std::atomic<uint32_t>> auxValueVector(nclus); // auxiliary vector that will contain all the non-zero values of each cluster (element of rowM)

    std::vector<reductionThread> redArgs;
    redArgs.reserve(numThreadsClus);

    std::vector<sumThread> sumArgs;
    sumArgs.reserve(numThreads);

    std::vector<assignThread> assignArgs;
    assignArgs.reserve(numThreadsClus);

    for (size_t id = 1; id < (nclus + 1); id++) // cluster ids start from 1
    {
        rowM[id - 1] = allCount;
        clusterHasElements = 0;

        auxValueVector.assign(nclus * numThreads, 0);
        sumArgs.clear();

        for (size_t i = 0; i < numThreads; i++)
        {
            size_t start = i * chunk;
            size_t end = (i == (numThreads - 1)) ? (start + lastThreadChunk) : (start + chunk);
            sumArgs.push_back({id, start, end, row, col, val, c, auxValueVector, i * nclus, clusterHasElements});
            pthread_create(&threadsSum[i], NULL, fnSumAux, (void *)&sumArgs[i]);
        }

        for (pthread_t i : threadsSum)
        {
            pthread_join(i, NULL);
        }

        if (!clusterHasElements)
            continue;

        for (size_t i = 0; i < numThreadsClus; i++)
        {
            size_t id_start = i * chunkClus + 1;
            size_t id_end = (i == (numThreadsClus - 1)) ? (id_start + lastThreadChunkClus) : (id_start + chunkClus);
            // printf("id_start, id_end: %ld, %ld\n", id_start, id_end);
            redArgs.push_back({id_start, id_end, nclus, auxValueVector, valMVector});
            pthread_create(&reductionThreads[i], NULL, fnReduction, (void *)&redArgs[i]);
        }

        for (pthread_t i : reductionThreads)
        {
            pthread_join(i, NULL);
        }

        // for (size_t i = 0; i < nclus * numThreads; i++)
        // {
        //     valMVector[(i % nclus)] += auxValueVector[i]; // reduction
        // }

        for (size_t i = 0; i < numThreadsClus; i++)
        {
            size_t start = i * chunkClus;
            size_t end = (i == (numThreadsClus - 1)) ? (start + lastThreadChunkClus) : (start + chunkClus);
            assignArgs.push_back({start, end, valMVector, allCount, colM, valM});
            pthread_create(&threadsAssign[i], NULL, fnAssignM, (void *)&assignArgs[i]);
        }

        for (pthread_t i : threadsAssign)
        {
            pthread_join(i, NULL);
        }

        assignArgs.clear();
        redArgs.clear();
        valMVector.assign(nclus, 0);
    }

    colM.resize(allCount.load());
    valM.resize(allCount.load());
}