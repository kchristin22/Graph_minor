#pragma once

#include <omp.h>
#include "coo_to_csr.hpp"

#define ELEMENTS_PER_CACHE_LINE_INT (64 / sizeof(int))
#define ELEMENTS_PER_CACHE_LINE_SIZE_T (64 / sizeof(size_t))

inline void calChunk(size_t &chunk, const size_t n, const size_t min_chunk, const uint32_t numThreads)
{
    size_t cacheLines = n / min_chunk; // a chunk of x elements of a size_t array (8*x bytes) will be equal to a cache line (64 bytes)
    if (!cacheLines)                   // the array of the smallest type (int < size_t) fits in a cache line
    {
        cacheLines = 1;
    }

    chunk = cacheLines * min_chunk / numThreads; // min chunk size
    if (!chunk)                                  // if we have too many threads, then threads must share a cache line
    {
        chunk = n / numThreads; // we assign equal number of elements to each thread
    }
}

inline void numClusters(size_t &nclus, const std::vector<size_t> &c, const uint32_t numThreads);

void GMopenMP(CSR &csrM, const CSR &csr, const std::vector<size_t> &c, const uint32_t numThreads);