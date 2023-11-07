#include "sequential.hpp"
#include "stdio.h"
#include "cstdlib"
#include "cmath"

void seq(std::vector<int> &M, size_t *L, std::vector<int> &A, std::vector<size_t> &c)
{
    printf("Hello from seq\n");

    if (A.size() != (c.size() * c.size()))
    {
        printf("Error: sizes of A and c are incompatible\n");
        exit(1);
    }
    size_t n = c.size();

    size_t n2 = (size_t)pow(n, 2);
    std::vector<size_t> row(n);             // offset of val vector for each row
    std::vector<size_t> col(n2);            // cluster id of the column of each element of val vector
    std::vector<int> val(n2);               // vector containing the non-zero elements (negative of positive) of A
    std::vector<size_t> discreetClus(n, 0); // vector where the ith element is a if cluster i has a nodes

    size_t allCount = 0, x; // volatile in parallel maybe?

    for (size_t i = 0; i < n; i++)
    {
        discreetClus[c[i]] += 1; // we assume that there is no ith row and column that are both zero so we know that all ids included in c exist in A
                                 // by atomically adding 1 to the cluster of the ith row we know how many nodes are in each cluster

        row[i] = allCount; // store offset of val vector for each row
        for (size_t j = 0; j < n; j++)
        {
            x = A[i * n + j];
            if (x == 0)
                continue;

            val[allCount] = x;
            col[allCount] = c[j]; // check if ids start from 1 and are continuous
            allCount++;
        }
    }

    val.resize(allCount); // resize vectors to the actual size of the non-zero elements to be more cache-efficient
    col.resize(allCount);

    size_t nclus = 0;
    for (size_t i = 0; i < n; i++) // check if ids start from 1 and are continuous
    {
        if (discreetClus[i] == 0)
            continue;
        nclus += 1;
    }
    *L = nclus;

    size_t end;
    std::vector<int> colCompressed(n * nclus, 0);

    for (size_t i = 0; i < n; i++)
    {
        if (i == (n - 1))
            end = allCount;
        else
            end = row[i + 1];
        for (size_t j = row[i]; j < end; j++)
        {
            colCompressed[i * nclus + col[j]] += val[j]; // compress cols by summing the values of each cluster to the column the cluster id points to
        }
    }

    // if ids are continuous they don't need sorting

    M.resize(nclus * nclus);

    for (size_t id = 0; id < nclus; id++) // check if ids start from 1 and are continuous
    {
        for (size_t i = 0; i < n; i++)
        {
            if (id != c[i]) // c[i]: cluster of row i of colCompressed/row
                continue;

            for (size_t j = 0; j < nclus; j++)
            {
                M[id * nclus + j] += colCompressed[i * nclus + j]; // compress rows of colCompressed by summing the rows of each cluster to the row of M the cluster id points to
            }
        }
    }
};