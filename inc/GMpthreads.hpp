#pragma once

#include <pthread.h>
#include <atomic>
#include <queue>
#include "coo_to_csr.hpp"

#define ELEMENTS_PER_CACHE_LINE_INT (64 / sizeof(int)) // size of a cache line is 64 bytes
#define ELEMENTS_PER_CACHE_LINE_SIZE_T (64 / sizeof(size_t))

/* Arguments of the threads responsible to calculate the number of clusters */
struct nclusThread
{
    const size_t start; // start of range of elements of a vector 1xn to be processed by the thread
    const size_t end;
    const std::vector<size_t> &c;
    size_t &nclus; // local max id measured by the thread
};

/* Arguments of the threads executing the main algorithm */
struct threadArgs
{
    const size_t nclus; // the overall number of clusters
    const size_t start; // start of range of elements of a vector 1xn to be processed by the thread
    const size_t end;
    const size_t startClus; // start of range of elements of a vector 1xnclus to be processed by the thread
    const size_t endClus;
    const CSR &csr; // the input matrix in CSR format
    const std::vector<size_t> &c;
    std::vector<std::atomic<uint32_t>> &commonAux; // common auxiliary vector to save compressed row and col of a single cluster
    std::atomic<size_t> &allCount;                 // common counter to keep track of the index of the col and val CSR output
    CSR &csrM;                                     // output matrix in CSR format
};

/* Structure only used in the pthreads version with dynamic scheduling. This is the task saved in the queue by the prodcuers and executed by consumers. */
struct Task
{
    std::queue<Task> *queue;
    pthread_mutex_t *queueMutex;
    pthread_cond_t *queueNotEmpty;
    size_t chunkStart;
    size_t chunkEnd;
    // void *args;
    nclusThread nclusArgs;
    void *(*fn)(void *);
};

/* Structure only used in the pthreads version with dynamic scheduling. This is the input arguments of a consumer thread.  */
struct consThread
{
    std::queue<Task> *queue;
    pthread_mutex_t *queueMutex;
    pthread_cond_t *queueNotEmpty;
};

/* Calculate the number of elements (chunk) to appoint to a thread to make computations thread-safe
 * @params
 *   chunk (output): the number of elements to appoint to a thread
 *   lastThreadChunk (output): the number of elements to appoint to the last thread, in case the division is not exact
 *   n (input): the number of elements to be divided into chunks
 *   min_chunk (input): the minimum number of elements to appoint to a thread
 *   numThreads (input): the number of threads to be used
 */
inline void calChunk(size_t &chunk, size_t &lastThreadChunk, const size_t n, const size_t min_chunk, const uint32_t numThreads);

/* Calculates the number of discrete clusters, in the case of having contiguous ids */
void *fnNumClusters(void *args);

/* The sequential algorithm implemented using pthreads */
void *fnThread(void *args);

/* Function that handles the threads used for computing the graph minor of the input graph
 * @params:
 *   csrM (output): CSR representation of the adjacency matrix corresponding to the graph minor of the input graph
 *   csr (input): CSR representation of the adjacency matrix of the input graph
 *   c (input): vector containing the cluster of each node of the input graph
 *   numThreads (input): the number of threads to be used
 */
void GMpthreads(CSR &csrM, const CSR &csr, const std::vector<size_t> &c, const uint32_t numThreads);
