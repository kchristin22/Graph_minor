#include "pthread.h"
#include "vector"
#include "stdio.h"
#include "atomic"
#include "stdint.h"
#include "queue"

#define ELEMENTS_PER_CACHE_LINE (64 / sizeof(int))

struct range
{
    size_t start;
    size_t end;
};

struct nclusThread
{
    size_t start;
    size_t end;
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

struct sumThread
{
    const size_t id;
    const size_t start;
    const size_t end;
    const std::vector<size_t> &row;
    const std::vector<size_t> &col;
    const std::vector<uint32_t> &val;
    const std::vector<size_t> &c;
    std::vector<uint32_t> &auxValueVector;
size_t auxValueIndex;
    bool &clusterHasElements;
};

struct assignThread
{
    const size_t start;
    const size_t end;
    std::vector<uint32_t> &auxValueVector;
    std::atomic<size_t> &allCount;
    std::vector<size_t> &colM;
    std::vector<uint32_t> &valM;
};

void *fnNumClusters(void *args);

void *fnClearAux(void *args);

void *fnSumAux(void *args);

void *fnAssignM(void *args);

void pthreads(std::vector<size_t> &rowM, std::vector<size_t> &colM, std::vector<uint32_t> &valM,
              const std::vector<size_t> &row, const std::vector<size_t> &col, const std::vector<uint32_t> &val, const std::vector<size_t> &c);