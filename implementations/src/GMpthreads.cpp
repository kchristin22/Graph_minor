#include <iostream>
#include "GMpthreads.hpp"

pthread_barrier_t barrier;

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

void *fnThread(void *args) // need numThreads, numClus, CSR, CSRM, id of thread (to calculate range), chunk, chunkClus, allCount, current clus id
{

    threadArgs *Args = (threadArgs *)args;
    size_t start = Args->start;
    size_t end = Args->end; // add cluster start index

    if (end > (Args->c.size()))
        end = Args->c.size();
    size_t startClus = Args->startClus;
    size_t endClus = Args->endClus;
    if (endClus > Args->nclus)
        endClus = Args->nclus;

    std::vector<uint32_t> auxValueVector(Args->nclus, 0);

    size_t x, endAux, localCount = 0;

    for (size_t id = 1; id < (Args->nclus + 1); id++)
    { // cluster ids start from 1

        pthread_barrier_wait(&barrier);

        Args->csrM.row[id - 1] = Args->allCount.load();

        auxValueVector.assign(Args->nclus, 0);

        for (size_t i = start; i < end; i++)
        {
            if (id != Args->c[i]) // c[i] : cluster of row i of colCompressed / row
                continue;

            endAux = Args->csr.row[i + 1];

            x = Args->csr.row[i];
            if (x == endAux) // this row has no non-zero elements
                continue;

            for (size_t j = x; j < endAux; j++)
            {
                auxValueVector[Args->c[Args->csr.col[j]] - 1] += (Args->csr.val[j]); // compress cols by summing the values of each cluster to the column the cluster id points to
            }
        }

        for (size_t i = 0; i < Args->nclus; i++)
        {
            if (auxValueVector[i] == 0)
                continue;
            Args->commonAux[i].fetch_add(auxValueVector[i]);
        }

        pthread_barrier_wait(&barrier);

        if (!Args->doAssign)
            continue;

        for (size_t i = startClus; i < endClus; i++)
        {
            if (Args->commonAux[i].load() == 0)
                continue;

            localCount = Args->allCount.fetch_add(1); // returns the previous value

            Args->csrM.val[localCount] = Args->commonAux[i];
            Args->csrM.col[localCount] = i;
            Args->commonAux[i] = 0;
        }
    }

    return nullptr;
}

void GMpthreads(CSR &csrM, const CSR &csr, const std::vector<size_t> &c)
{
    if (csr.row.size() != (c.size() + 1))
    {
        printf("Error: sizes of row and c are incompatible\n");
        exit(1);
    }
    else if (csr.col.size() != csr.val.size())
    {
        printf("Error: sizes of col and val are incompatible\n");
        exit(1);
    }
    else if (csr.row.size() == (csr.col.size() + 1))
    {
        printf("Error: CSR requires more space than dense matrix representation \n Use dense matrix implementation instead...\n");
        exit(1);
    }
    else if ((csrM.row.size() != csr.row.size()) || (csrM.col.size() != csr.col.size()) || (csrM.val.size() != csr.val.size()))
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

    if (n > (chunk * numThreads)) // n does not fill in an integer number of cachelines
    {
        lastThreadChunk = chunk + (n - (chunk * numThreads)); // the last thread gets its chunk + the remaining items
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

    for (pthread_t &i : nclusThreads)
    {
        pthread_join(i, NULL);
    }

    for (size_t i = 0; i < numThreads; i++)
    {
        nclus += nclusVector[i];
    }

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
    // numThreadsClus = 4;

    size_t chunkClus = cacheLinesClus * ELEMENTS_PER_CACHE_LINE / numThreadsClus;
    size_t lastThreadChunkClus = chunkClus;

    if (nclus > (chunkClus * numThreadsClus)) // n does not fill in an integer number of cachelines
    {
        lastThreadChunkClus = chunkClus + (nclus - (chunkClus * numThreadsClus)); // the last thread gets its chunk + the remaining items
    }
    else if (nclus <= (cacheLinesClus * ELEMENTS_PER_CACHE_LINE)) // the array fits in a single cache line
        lastThreadChunkClus = nclus;

    std::atomic<size_t> allCount(0);
    csrM.row.resize(nclus + 1); // resize vector to the number of clusters

    std::cout << "Running with numThreads, numThreadsClus: " << numThreads << ", " << numThreadsClus << std::endl;

    pthread_t threads[numThreads];

    std::vector<threadArgs> argsVector;
    argsVector.reserve(numThreads);

    bool doAssign = 1;

    std::vector<std::atomic<uint32_t>> commonAux(nclus); // auxiliary vector that will contain all the non-zero values of each cluster (element of rowM)
    pthread_barrier_init(&barrier, NULL, numThreads);

    size_t start, startClus;

    for (size_t i = 0; i < numThreads; i++)
    {
        start = i * chunk;
        size_t end = (i == (numThreads - 1)) ? (start + lastThreadChunk) : (start + chunk);
        startClus = i * chunkClus;
        size_t endClus = (i == (numThreadsClus - 1)) ? (startClus + lastThreadChunkClus) : (startClus + chunkClus);
        if (i == numThreadsClus)
            doAssign = 0;

        argsVector.push_back({nclus, start, end, startClus, endClus, doAssign, csr, c, commonAux, allCount, csrM});
        pthread_create(&threads[i], NULL, fnThread, (void *)&argsVector[i]);
    }

    for (pthread_t &i : threads)
    {
        pthread_join(i, NULL);
    }

    csrM.row[nclus] = allCount; // last element of rowM is the number of non-zero elements of valM and colM
    csrM.col.resize(allCount);
    csrM.val.resize(allCount);

    pthread_barrier_destroy(&barrier);
}