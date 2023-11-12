#include "openMP.hpp"
#include "cstdlib"

inline void numClusters(size_t &nclus, std::vector<size_t> &c)
{
    size_t n = c.size();
    std::vector<size_t> discreetClus(n, 0); // vector where the ith element is a if cluster i has a nodes

    size_t chunk = n / 8; // a chunk of x elements of a size_t array (8*x bytes) will be equal to a cache line (64 bytes)
    size_t numThreads = (chunk / 4);
    if (!chunk)
    {
        chunk = n;
    }

    if (!numThreads)
    {
        numThreads = 1;
    }
    else if (numThreads > 4)
    {
        numThreads = 4;
    }

#pragma omp parallel num_threads(numThreads)
#pragma omp for nowait schedule(static, chunk)
    for (size_t i = 0; i < n; i++)
    {
        discreetClus[c[i] - 1] = 1; // we assume that there is no ith row and column that are both zero so we know that all ids included in c exist in A
                                    // we can atomically add 1, instead, to the cluster of the ith row to know how many nodes are in each cluster
    }

    nclus = 0;

#pragma omp parallel num_threads(numThreads)
#pragma omp for nowait reduction(+ : nclus) schedule(dynamic, chunk)
    for (size_t i = 0; i < n; i++)
    {
        if (discreetClus[i] == 0)
            continue;
        nclus += 1;
    }
}

void openMP(std::vector<size_t> &rowM, std::vector<size_t> &colM, std::vector<int> &valM,
            std::vector<size_t> &row, std::vector<size_t> &col, std::vector<int> &val, std::vector<size_t> &c)
{
    if (row.size() != c.size())
    {
        printf("Error: sizes of row and c are incompatible\n");
        exit(1);
    }
    else if (col.size() != val.size())
    {
        printf("Error: sizes of col and val are incompatible\n");
        exit(1);
    }
    else if (row.size() == col.size())
    {
        printf("Error: CSR requires more space than dense matrix representation \n Use dense matrix implementation instead...\n");
        exit(1);
    }
    else if ((rowM.size() != row.size()) || (colM.size() != col.size()) || (valM.size() != val.size()))
    // if the input graph is its the minor, the dimensions of the compressed vectors should be equal to the dimensions of the input vectors
    {
        printf("Error: at least one of the compressed vectors doesn't have enough allocated space \n");
        exit(1);
    }

    size_t n = c.size();
    size_t nz = val.size();
    size_t nclus;

    numClusters(nclus, c); // find the number of distinct clusters

    size_t end, allCount = 0, localCount;
    bool clusterHasElements = 0;
    int auxValueVector[nclus]{0}; // auxiliary vector that will contain all the non-zero values of each cluster (element of rowM)
    rowM.resize(nclus);           // resize vector to the number of clusters

    size_t chunk = n / 16; // how many cache lines the array fills:
                           // a chunk of x elements of an int array (4*x bytes) will be equal to a cache line (64 bytes) to avoid false sharing
    size_t numThreads = chunk / 4;
    if (!chunk) // the array of the smallest type (int < size_t) fits in a cache line
    {
        chunk = n;
        numThreads = 1;
    }

    if (!numThreads)
    {
        numThreads = 1;
    }
    else if (numThreads > 4)
    {
        numThreads = 4;
    }

    size_t chunkClus = nclus / 16;
    size_t numThreadsClus = chunkClus / 4;
    if (!chunkClus) // the auxValueVector array fits in a cache line
    {
        chunkClus = nclus;
    }

    if (!numThreadsClus)
    {
        numThreadsClus = 1;
    }
    else if (numThreadsClus > 4)
    {
        numThreadsClus = 4;
    }

    for (size_t id = 1; id < (nclus + 1); id++) // cluster ids start from 1
    {
        rowM[id - 1] = allCount;

#pragma omp parallel num_threads(numThreadsClus)
#pragma omp for nowait schedule(static, chunkClus)
        for (size_t i = 0; i < nclus; i++)
            auxValueVector[i] = 0; // reset auxiliary vector

        clusterHasElements = 0;

#pragma omp parallel num_threads(numThreads)
#pragma omp for nowait reduction(+ : auxValueVector[ : nclus]) private(end) schedule(dynamic, chunk)
        // reduction of each element of the auxiliary vector
        for (size_t i = 0; i < n; i++)
        {
            if (id != c[i]) // c[i]: cluster of row i of colCompressed/row
                continue;

            if (i == (n - 1)) // last row contains the last range of non-zero elements
                end = nz;
            else
                end = row[i + 1];

            if (row[i] == end) // this row has no non-zero elements
                continue;
            if (!clusterHasElements) // declare that this row has non-zero elements only once
                clusterHasElements = 1;

            for (size_t j = row[i]; j < end; j++)
            {
                auxValueVector[c[col[j]] - 1] += val[j]; // compress cols by summing the values of each cluster to the column the cluster id points to
            }
        }

        if (!clusterHasElements)
            continue;

#pragma omp parallel num_threads(numThreadsClus)
#pragma omp for nowait private(localCount) schedule(dynamic, chunkClus)
        for (size_t i = 0; i < nclus; i++)
        {
            if (auxValueVector[i] == 0)
                continue;
#pragma omp atomic capture
            {
                localCount = allCount;
                allCount++;
            }
            valM[localCount] = auxValueVector[i];
            colM[localCount] = i;
        }
    }

    colM.resize(allCount);
    valM.resize(allCount);
}