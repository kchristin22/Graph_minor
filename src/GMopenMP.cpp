#include "GMopenMP.hpp"
#include <pthread.h>

inline void numClusters(size_t &nclus, const std::vector<size_t> &c, const uint32_t numThreads)
{
    size_t n = c.size();

    size_t chunk;
    calChunk(chunk, n, ELEMENTS_PER_CACHE_LINE_SIZE_T, numThreads);

    size_t maxVal = 0;

#pragma omp parallel num_threads(numThreads)
#pragma omp for nowait schedule(static, chunk) reduction(max : maxVal)
    for (size_t i = 0; i < n; i++)
    {
        maxVal = (c[i] > maxVal) ? c[i] : maxVal;
    }
    nclus = maxVal;
}

void GMopenMP(CSR &csrM, const CSR &csr, const std::vector<size_t> &c, const uint32_t numThreads)
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
    size_t nclus;

    numClusters(nclus, c, numThreads); // find the number of distinct clusters

    size_t end, allCount = 0, localCount;
    bool clusterHasElements = 0;
    uint32_t auxValueVector[nclus]{0}; // auxiliary vector that will contain all the non-zero values of each cluster (element of rowM)
    csrM.row.resize(nclus + 1);        // resize vector to the number of clusters

    size_t chunk, chunkClus;
    calChunk(chunk, n, ELEMENTS_PER_CACHE_LINE_INT, numThreads); // how many cache lines the array fills:
                                                                 // a chunk of x elements of an int array (4*x bytes) will be equal to a cache line (64 bytes) to avoid false sharing
                                                                 // the array of the smallest type (int < size_t) fits in a cache line

    calChunk(chunkClus, nclus, ELEMENTS_PER_CACHE_LINE_INT, numThreads); // how many cache lines the array fills:
                                                                         // a chunk of x elements of an int array (4*x bytes) will be equal to a cache line (64 bytes) to avoid false sharing
                                                                         // the array of the smallest type (int < size_t) fits in a cache line

    for (size_t id = 1; id < (nclus + 1); id++) // cluster ids start from 1
    {
        csrM.row[id - 1] = allCount;

#pragma omp parallel num_threads(numThreads)
#pragma omp for nowait schedule(static, chunkClus)
        for (size_t i = 0; i < nclus; i++)
        {
            auxValueVector[i] = 0; // reset auxiliary vector
            // printf("thread id: %d ", omp_get_thread_num());
        }

        // printf("thread id: %ld\n", pthread_self());
        // printf("exited, main thread id: %d ", omp_get_thread_num());

        clusterHasElements = 0;

#pragma omp parallel num_threads(numThreads)
#pragma omp for nowait reduction(+ : auxValueVector[ : nclus]) private(end) schedule(dynamic, chunk)
        // reduction of each element of the auxiliary vector
        for (size_t i = 0; i < n; i++)
        {
            // printf("thread id: %d ", omp_get_thread_num());
            if (id != c[i]) // c[i]: cluster of row i of colCompressed/row
                continue;

            end = csr.row[i + 1];

            if (csr.row[i] == end) // this row has no non-zero elements
                continue;
            if (!clusterHasElements) // declare that this row has non-zero elements only once
                clusterHasElements = 1;

            for (size_t j = csr.row[i]; j < end; j++)
            {
                auxValueVector[c[csr.col[j]] - 1] += csr.val[j]; // compress cols by summing the values of each cluster to the column the cluster id points to
            }
        }

        if (!clusterHasElements)
            continue;

#pragma omp parallel num_threads(numThreads)
#pragma omp for nowait private(localCount) schedule(dynamic, chunkClus)
        for (size_t i = 0; i < nclus; i++)
        {
            if (auxValueVector[i] == 0)
                continue;
#pragma omp atomic capture
            {
                localCount = allCount;
                allCount++;
            }
            csrM.val[localCount] = auxValueVector[i];
            csrM.col[localCount] = i;
        }
    }

    csrM.row[nclus] = allCount; // last element of rowM is the number of non-zero elements of valM and colM
    csrM.col.resize(allCount);
    csrM.val.resize(allCount);
}