#include "pthread.h"
#include "vector"
#include "stdio.h"
#include "atomic"

struct forThread
{
    size_t start;
    size_t end;
    std::vector<std::atomic<int>> &array;
};

struct sumThread
{
    size_t id;
    size_t start;
    size_t end;
    std::vector<size_t> &row;
    std::vector<size_t> &col;
    std::vector<int> &val;
    std::vector<size_t> &c;
    std::vector<std::atomic<int>> &auxValueVector;
    bool &clusterHasElements;
};

struct assignThread
{
    size_t start;
    size_t end;
    std::vector<std::atomic<int>> &auxValueVector;
    std::atomic<size_t> &allCount;
    std::vector<size_t> &colM;
    std::vector<int> &valM;
};

inline void numClusters(size_t &nclus, std::vector<size_t> &c);

inline void *fnClearAux(void *args);

void *fnSumAux(void *args);

void *fnAssignM(void *args);

void pthreads(std::vector<size_t> &rowM, std::vector<size_t> &colM, std::vector<int> &valM,
              std::vector<size_t> &row, std::vector<size_t> &col, std::vector<int> &val, std::vector<size_t> &c);