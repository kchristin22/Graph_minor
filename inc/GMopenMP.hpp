#pragma once

#include <omp.h>
#include "coo_to_csr.hpp"

#define ELEMENTS_PER_CACHE_LINE_INT (64 / sizeof(int)) // size of a cache line is 64 bytes
#define ELEMENTS_PER_CACHE_LINE_SIZE_T (64 / sizeof(size_t))

/* Calculate the number of elements (chunk) to appoint to a thread to make computations thread-safe. This function is used for coo-to-csr conversion as well.
 * @params
 *   chunk (output): the number of elements to appoint to a thread
 *   lastThreadChunk (output): the number of elements to appoint to the last thread, in case the division is not exact
 *   n (input): the number of elements to be divided into chunks
 *   min_chunk (input): the minimum number of elements to appoint to a thread
 *   numThreads (input): the number of threads
 */
inline void calChunk(size_t &chunk, const size_t n, const size_t min_chunk, const uint32_t numThreads)
{
    size_t cacheLines = n / min_chunk; // we assume that the minimum chunk given fills at least a whole cache line
                                       // f.i. a min chunk equal to 8 elements of a size_t array (8*8 bytes) will be equal to a cache line (64 bytes)

    if (!cacheLines) // the array fits in a single cache line
    {
        cacheLines = 1;
    }

    chunk = cacheLines * min_chunk / numThreads; // divide cache lines into chunks of equal size
    if (!chunk)                                  // if we have too many threads, then threads must share a cache line
    {
        chunk = n / numThreads; // we assign equal number of elements to each thread
    }

    // no need to calculate the remaining elements, OpenMP will take care of it
}

/* Calculates the number of discrete clusters. This function is an idea of what counting the ids would look like if they weren't contiguous.
 * @params:
 *   nclus (output): number of discrete clusters
 *   c (input): vector containing the cluster of each node of the input graph
 *   numThreads (input): the number of threads to be used
 */
inline void numClusters(size_t &nclus, const std::vector<size_t> &c, const uint32_t numThreads);

/* The sequential algorithm implemented using OpenMP
 * @params:
 *   csrM (output): CSR representation of the adjacency matrix corresponding to the graph minor of the input graph
 *   csr (input): CSR representation of the adjacency matrix of the input graph
 *   c (input): vector containing the cluster of each node of the input graph
 *   numThreads (input): the number of threads to be used
 */
void GMopenMP(CSR &csrM, const CSR &csr, const std::vector<size_t> &c, const uint32_t numThreads);