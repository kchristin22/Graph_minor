#pragma once

#include <pthread.h>
#include <atomic>
#include <queue>
#include "coo_to_csr.hpp"

#define ELEMENTS_PER_CACHE_LINE (64 / sizeof(int))

struct range
{
    size_t start;
    size_t end;
};

struct nclusThread
{
    const size_t start;
    const size_t end;
    const std::vector<size_t> &c;
    size_t &nclus;
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

struct forThread
{
    const size_t start;
    const size_t end;
    std::vector<std::atomic<uint32_t>> &array;
};

struct reductionThread
{
    const size_t id_start;
    const size_t id_end;
    const size_t nclus;
    const std::vector<uint32_t> &inputArray;
    std::vector<uint32_t> &outputArray;
};

struct sumThread
{
    const size_t id;
    const size_t start;
    const size_t end;
    const CSR &csr;
    const std::vector<size_t> &c;
    std::vector<uint32_t> &auxValueVector;
    size_t auxValueIndex;
    bool &clusterHasElements;
};

struct assignThread
{
    const size_t start;
    const size_t end;
    const std::vector<uint32_t> &auxValueVector;
    std::atomic<size_t> &allCount;
    CSR &csrM;
};

struct threadArgs
{
    const size_t threadID;
    const size_t nclus;
    const size_t start;
    const size_t chunk;
    const size_t startClus;
    const size_t chunkClus;
    const bool doAssign;
    const CSR &csr;
    const std::vector<size_t> &c;
    std::vector<std::atomic<uint32_t>> &commonAux;
    std::atomic<size_t> &allCount;
    CSR &csrM;
};

void *fnNumClusters(void *args);

void *fnClearAux(void *args);

void *fnSumAux(void *args);

void *fnAssignM(void *args);

void GMpthreads(CSR &csrM, const CSR &csr, const std::vector<size_t> &c);
