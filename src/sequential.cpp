#include "sequential.hpp"
#include "stdio.h"
#include "cstdlib"

inline void numClusters(size_t &nclus, std::vector<size_t> &c)
{
    size_t n = c.size();
    std::vector<size_t> discreetClus(n, 0); // vector where the ith element is a if cluster i has a nodes

    for (size_t i = 0; i < n; i++)
    {
        discreetClus[c[i] - 1] = 1; // we assume that there is no ith row and column that are both zero so we know that all ids included in c exist in A
                                    // we can atomically add 1, instead, to the cluster of the ith row to know how many nodes are in each cluster
    }

    nclus = 0;
    for (size_t i = 0; i < n; i++)
    {
        if (discreetClus[i] == 0)
            continue;
        nclus += 1;
    }
}

void seq(std::vector<int> &M, std::vector<int> &A, std::vector<size_t> &c)
{
    printf("Hello from seq\n");

    if (A.size() != (c.size() * c.size()))
    {
        printf("Error: sizes of A and c are incompatible\n");
        exit(1);
    }
    size_t n = c.size();
    size_t nclus;

    numClusters(nclus, c);

    int val;
    std::vector<int> colCompressed(n * nclus, 0);

    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            val = A[i * n + j];
            if (val == 0)
                continue;
            colCompressed[i * nclus + c[j] - 1] += val; // compress cols by summing the values of each cluster to the column the cluster id points to
        }
    }

    // if ids are continuous they don't need sorting

    M.resize(nclus * nclus, 0);

    for (size_t id = 1; id < (nclus + 1); id++) // cluster ids start from 1
    {                                           // access one row of M at a time
        for (size_t i = 0; i < n; i++)
        {
            if (id != c[i]) // c[i]: cluster of row i of colCompressed/row
                continue;

            for (size_t j = 0; j < nclus; j++)
            {
                M[(id - 1) * nclus + j] += colCompressed[i * nclus + j]; // compress rows of colCompressed by summing the rows of each cluster to the row of M the cluster id points to
            }
        }
    }
};

void seq(std::vector<size_t> &rowM, std::vector<size_t> &colM, std::vector<int> &valM,
         std::vector<size_t> &row, std::vector<size_t> &col, std::vector<int> &val, std::vector<size_t> &c)
{
    printf("Hello from seq\n");

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

    numClusters(nclus, c);

    size_t end, clusterHasElements, allCount = 0; // store offset to assign to each rowM element
    std::vector<int> auxValueVector(nclus);       // auxiliary vector that will contain all the non-zero values of each cluster (element of rowM)
    rowM.resize(nclus);                           // resize vector to the number of clusters

    for (size_t id = 1; id < (nclus + 1); id++) // cluster ids start from 1
    {
        rowM[id - 1] = allCount;
        auxValueVector.assign(nclus, 0); // reset auxiliary vector
        clusterHasElements = 0;

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
        for (size_t i = 0; i < nclus; i++)
        {
            if (auxValueVector[i] == 0)
                continue;

            valM[allCount] = auxValueVector[i];
            colM[allCount] = i;
            allCount++;
        }
    }
    colM.resize(allCount);
    valM.resize(allCount);
}

void seqDenseCSR(std::vector<size_t> &rowM, std::vector<size_t> &colM, std::vector<int> &valM, std::vector<size_t> &row, std::vector<size_t> &col, std::vector<int> &val, std::vector<size_t> &c)
{
    printf("Hello from seq\n");

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

    numClusters(nclus, c);

    size_t end;
    std::vector<int> colCompressed(n * nclus, 0);

    for (size_t i = 0; i < n; i++)
    {
        if (i == (n - 1))
            end = nz;
        else
            end = row[i + 1];
        for (size_t j = row[i]; j < end; j++)
        {
            colCompressed[i * nclus + c[col[j]] - 1] += val[j]; // compress cols by summing the values of each cluster to the column the cluster id points to
        }
    }

    // if ids are continuous they don't need sorting

    std::vector<int> M(nclus * nclus, 0);

    printf("Here \n");

    for (size_t id = 1; id < (nclus + 1); id++) // cluster ids start from 1
    {                                           // access one row of M at a time
        for (size_t i = 0; i < n; i++)
        {
            if (id != c[i]) // c[i]: cluster of row i of colCompressed/row
                continue;

            for (size_t j = 0; j < nclus; j++)
            {
                M[(id - 1) * nclus + j] += colCompressed[i * nclus + j]; // compress rows of colCompressed by summing the rows of each cluster to the row of M the cluster id points to
            }
        }
    }

    size_t count = 0;
    int x;
    rowM.resize(nclus);
    for (size_t i = 0; i < nclus; i++)
    {
        rowM[i] = count;
        for (size_t j = 0; j < nclus; j++)
        {
            x = M[i * nclus + j];
            if (x == 0)
                continue;
            valM[count] = x;
            colM[count] = j;
            count++;
        }
    }
    colM.resize(count);
    valM.resize(count);
}