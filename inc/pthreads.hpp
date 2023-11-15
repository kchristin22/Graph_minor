#include "pthread.h"
#include "vector"
#include "stdio.h"
#include "atomic"
#include "stdint.h"

#define ELEMENTS_PER_CACHE_LINE (64 / sizeof(int))

struct nclusThread
{
    const size_t start;
    const size_t end;
    const std::vector<size_t> &c;
    std::atomic<size_t> &nclus;
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
    std::vector<std::atomic<uint32_t>> &auxValueVector;
    bool &clusterHasElements;
};

struct assignThread
{
    const size_t start;
    const size_t end;
    std::vector<std::atomic<uint32_t>> &auxValueVector;
    std::atomic<size_t> *allCount;
    std::vector<size_t> &colM;
    std::vector<uint32_t> &valM;
};

inline void *fnNumClusters(void *args);

inline void *fnClearAux(void *args);

void *fnSumAux(void *args);

void *fnAssignM(void *args);

void pthreads(std::vector<size_t> &rowM, std::vector<size_t> &colM, std::vector<uint32_t> &valM,
              const std::vector<size_t> &row, const std::vector<size_t> &col, const std::vector<uint32_t> &val, const std::vector<size_t> &c);