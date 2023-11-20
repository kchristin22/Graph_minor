#pragma once

#include <pthread.h>
#include <atomic>
#include <queue>
#include "coo_to_csr.hpp"

#define ELEMENTS_PER_CACHE_LINE_INT (64 / sizeof(int))
#define ELEMENTS_PER_CACHE_LINE_SIZE_T (64 / sizeof(size_t))

struct nclusThread
{
    const size_t start;
    const size_t end;
    const std::vector<size_t> &c;
    size_t &nclus;
};

struct threadArgs
{
    const size_t nclus;
    const size_t start;
    const size_t end;
    const size_t startClus;
    const size_t endClus;
    const CSR &csr;
    const std::vector<size_t> &c;
    std::vector<std::atomic<uint32_t>> &commonAux;
    std::atomic<size_t> &allCount;
    CSR &csrM;
};

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

struct consThread
{
    std::queue<Task> *queue;
    pthread_mutex_t *queueMutex;
    pthread_cond_t *queueNotEmpty;
};

inline void calChunk(size_t &chunk, size_t &lastThreadChunk, const size_t n, const size_t min_chunk, const uint32_t numThreads);

void *fnNumClusters(void *args);

void GMpthreads(CSR &csrM, const CSR &csr, const std::vector<size_t> &c, const uint32_t numThreads);
