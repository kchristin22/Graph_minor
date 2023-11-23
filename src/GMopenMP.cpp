#include "GMopenMP.hpp"
#include <pthread.h>

inline void numClusters(size_t &nclus, const std::vector<size_t> &c, const uint32_t numThreads)
{
    size_t n = c.size();

    size_t chunk;
    calChunk(chunk, n, ELEMENTS_PER_CACHE_LINE_SIZE_T, numThreads); // c vector elements are size_t -> min chunk = cache line / sizeof(size_t)

    size_t maxVal = 0;

// perform a max reduction for cluster ids
#pragma omp parallel num_threads(numThreads)
#pragma omp for nowait schedule(static, chunk) reduction(max : maxVal) // static scheduling as we expect the workload to be equal for each thread
                                                                       // nowait allows threads to suspend after finishing their work in the loop,
                                                                       // without having to wait for the others to finish
    for (size_t i = 0; i < n; i++)
    {
        maxVal = (c[i] > maxVal) ? c[i] : maxVal;
    }
    nclus = maxVal;
}

void GMopenMP(CSR &csrM, const CSR &csr, const std::vector<size_t> &c, const uint32_t numThreads)
{
    size_t n = c.size();
    size_t nz = csr.val.size();

    // the row vector should be of size equal to the rows of the matrix, plus one (the latter is an implementation specific)
    // we assume that the c vector's size is set correctly equal to the number of nodes of the input grapg / rows of the input matrix
    if (csr.row.size() != (n + 1))
    {
        printf("Error: sizes of row and c are incompatible\n");
        exit(1);
    }
    else if (csr.col.size() != nz)
    {
        printf("Error: sizes of col and val are incompatible\n");
        exit(1);
    }
    else if (((n * n) - 2 * n) < (2 * nz)) // n^2 > n + 2*nz,  resulting in CSR taking more space than dense matrix representation
    {
        printf("Error: CSR requires more space than dense matrix representation \n Use dense matrix implementation instead...\n");
        exit(1);
    }
    else if ((csrM.row.size() != (n + 1)) || (csrM.col.size() != nz) || (csrM.val.size() != nz))
    // in case the input graph is its minor and no further compression is possible, the dimensions of the compressed vectors should be equal to the dimensions of the input vectors
    {
        printf("Error: at least one of the compressed vectors doesn't have enough allocated space \n");
        exit(1);
    }

    size_t nclus;

    numClusters(nclus, c, numThreads); // find the number of distinct clusters

    size_t end, allCount = 0, localCount; // loop variables
    bool clusterHasElements = 0;          // flag that denotes if the cluster has any elements
    uint32_t auxValueVector[nclus]{0};    // auxiliary array that will contain all the non-zero values of each cluster (element of rowM)
    csrM.row.resize(nclus + 1);           // resize vector to the number of clusters

    size_t chunk, chunkClus;
    calChunk(chunk, n, ELEMENTS_PER_CACHE_LINE_INT, numThreads); // how many cache lines the array fills:
                                                                 // the array of the smallest type (uint32_t < size_t) in that block of code fits in a cache line

    calChunk(chunkClus, nclus, ELEMENTS_PER_CACHE_LINE_INT, numThreads);

    for (size_t id = 1; id < (nclus + 1); id++) // cluster ids start from 1
    {
        csrM.row[id - 1] = allCount;

#pragma omp parallel num_threads(numThreads)
#pragma omp for nowait schedule(static, chunkClus) // static scheduling as we expect the workload to be equal for each thread
                                                   // nowait allows threads to suspend after finishing their work in the loop,
                                                   // without having to wait for the others to finish
        for (size_t i = 0; i < nclus; i++)
        {
            auxValueVector[i] = 0; // reset auxiliary vector
        }

        clusterHasElements = 0;

// reduction of each element of the auxiliary vector
#pragma omp parallel num_threads(numThreads)
#pragma omp for nowait reduction(+ : auxValueVector[ : nclus]) private(end) schedule(dynamic, chunk)
        for (size_t i = 0; i < n; i++)
        {
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

        if (!clusterHasElements) // executed by the main thread
            continue;

#pragma omp parallel num_threads(numThreads)
#pragma omp for nowait private(localCount) schedule(dynamic, chunkClus)
        for (size_t i = 0; i < nclus; i++)
        {
            if (auxValueVector[i] == 0)
                continue;
#pragma omp atomic capture
            {
                localCount = allCount; // read the value of allCount and increment it atomically
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