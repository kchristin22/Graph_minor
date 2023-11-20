#include <atomic>
#include <cilk/opadd_reducer.h>
#include <cilk/cilk_api.h>
#include "GMopenCilk.hpp"

inline void zero_s(void *v) { *(size_t *)v = 0; }
inline void plus_s(void *l, void *r) { *(size_t *)l += *(size_t *)r; }

inline void numClusters(size_t cilk_reducer(zero_s, plus_s) & nclus, const std::vector<size_t> &c, const uint32_t numThreads)
{
    size_t n = c.size();
    std::vector<size_t> discreetClus(n, 0); // vector where the ith element is a if cluster i has a nodes

    size_t cacheLines = n / ELEMENTS_PER_CACHE_LINE_SIZE_T; // a chunk of x elements of a size_t array (8*x bytes) will be equal to a cache line (64 bytes)

    if (!cacheLines) // the array of the smallest type (int < size_t) fits in a cache line
    {
        cacheLines = 1;
    }

    size_t chunk = cacheLines * ELEMENTS_PER_CACHE_LINE_SIZE_T / numThreads; // min chunk size
    if (!chunk)                                                              // if we have too many threads, then threads must share a cache line
    {
        chunk = n / numThreads; // we assign equal number of elements to each thread
    }

#pragma cilk_grainsize = chunk
    cilk_for(size_t i = 0; i < n; i++)
    {
        // printf("id: %lu\n", __cilkrts_running_on_workers());
        discreetClus[c[i] - 1] = 1; // we assume that there is no ith row and column that are both zero so we know that all ids included in c exist in A
                                    // we can atomically add 1, instead, to the cluster of the ith row to know how many nodes are in each cluster
    }

    nclus = 0;

#pragma cilk_grainsize = chunk
    cilk_for(size_t i = 0; i < n; i++)
    {
        if (discreetClus[i] == 0) // benefit from predicting
            continue;
        nclus += 1;
    }
}

void GMopenCilk(CSR &csrM, const CSR &csr, const std::vector<size_t> &c, const uint32_t numThreads)
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
    size_t cilk_reducer(zero_s, plus_s) nclus = 0;

    numClusters(nclus, c, numThreads); // find the number of distinct clusters

    size_t end, localCount;
    std::atomic<uint32_t> allCount(0);
    bool clusterHasElements = 0;
    std::vector<std::atomic<uint32_t>> auxValueVector(nclus); // auxiliary vector that will contain all the non-zero values of each cluster (element of rowM)
    // std::vector<uint32_t> auxValueVector(nclus);

    csrM.row.resize(nclus + 1); // resize vector to the number of clusters

    size_t cacheLines = n / ELEMENTS_PER_CACHE_LINE_INT; // how many cache lines the array fills:
                                                         // a chunk of x elements of an int array (4*x bytes) will be equal to a cache line (64 bytes) to avoid false sharing
    if (!cacheLines)                                     // the array of the smallest type (int < size_t) fits in a cache line
    {
        cacheLines = 1;
    }

    size_t chunk = cacheLines * ELEMENTS_PER_CACHE_LINE_INT / numThreads; // min chunk size
    if (!chunk)                                                           // if we have too many threads, then threads must share a cache line
    {
        chunk = n / numThreads; // we assign equal number of elements to each thread
    }

    size_t cacheLinesClus = nclus / ELEMENTS_PER_CACHE_LINE_INT;
    if (!cacheLinesClus) // the auxValueVector array fits in a cache line
    {
        cacheLinesClus = 1;
    }

    size_t chunkClus = cacheLinesClus * ELEMENTS_PER_CACHE_LINE_INT / numThreads; // min chunk size
    if (!chunkClus)                                                               // if we have too many threads, then threads must share a cache line
    {
        chunkClus = nclus / numThreads; // we assign equal number of elements to each thread
    }

    for (size_t id = 1; id < (nclus + 1); id++) // cluster ids start from 1
    {
        csrM.row[id - 1] = allCount.load();

#pragma cilk_grainsize = cacheLinesClus
        cilk_for(size_t i = 0; i < nclus; i++)
            auxValueVector[i]
                .store(0); // reset auxiliary vector

        clusterHasElements = 0;


#pragma cilk_grainsize = cacheLines
        cilk_for(size_t i = 0; i < n; i++)
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
                auxValueVector[c[csr.col[j]] - 1] += (csr.val[j]); // compress cols by summing the values of each cluster to the column the cluster id points to
            }
        }

        if (!clusterHasElements)
            continue;

#pragma cilk_grainsize = cacheLinesClus
        cilk_for(size_t i = 0; i < nclus; i++)
        {
            if (auxValueVector[i] == 0)
                continue;

            localCount = allCount.fetch_add(1);
            csrM.val[localCount] = auxValueVector[i];
            csrM.col[localCount] = i;
        }
    }

    csrM.row[nclus] = allCount;
    csrM.col.resize(allCount);
    csrM.val.resize(allCount);
}