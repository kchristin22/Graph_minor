#pragma once

#include "iostream"
#include "stdio.h"
#include "cstdlib"
#include "vector"
#include "cmath"

template <size_t n>
void seq(int *M, size_t *L, size_t A[n][n], size_t c[n])
{
    printf("Hello from seq\n");

    // check in main:
    // if ((A[0].size() != n) || (A.size() != (size_t)(pow(n, 2))) || (c.size() != n))
    // {
    //     printf("Error: Incorrect size of A or c\n");
    //     exit(1);
    // }

    size_t n2 = (size_t)pow(n, 2);
    std::vector<size_t> row(n);             // offset of val vector for each row
    std::vector<size_t> col(n2);            // cluster id of the column of each element of val vector
    std::vector<int> val(n2);               // vector containing the non-zero elements (negative of positive) of A
    std::vector<size_t> discreetClus(n, 0); // vector where the ith element is a if cluster i has a nodes

    size_t allCount = 0, x; // volatile in parallel maybe?

    for (size_t i = 0; i < n; i++)
    {
        discreetClus[c[i] - 1] += 1; // we assume that there is no ith row and column that are both zero so we know that all ids included in c exist in A
                                     // by atomically adding 1 to the cluster of the ith row we know how many nodes are in each cluster

        row[i] = allCount; // store offset of val vector for each row
        for (size_t j = 0; j < n; j++)
        {
            x = A[i][j];
            if (x == 0)
                continue;
            val[allCount] = x;
            col[allCount] = c[j] - 1;
            allCount++;
        }
    }

    val.resize(allCount); // resize vectors to the actual size of the non-zero elements to be more cache-efficient
    col.resize(allCount);

    size_t nclus = 0;
    for (size_t i = 0; i < (n + 1); i++)
    {
        if (discreetClus[i] == 0)
            continue;
        nclus += 1;
    }
    *L = nclus;

    size_t end;
    int colCompressed[n][nclus] = {0};

    for (size_t i = 0; i < n; i++)
    {
        if (i == (n - 1))
            end = allCount;
        else
            end = row[i + 1];
        for (size_t j = row[i]; j < end; j++)
        {
            colCompressed[i][col[j]] += val[j]; // compress cols by summing the values of each cluster to the column the cluster id points to
        }
    }

    // if ids are continuous they don't need sorting

    M = (int *)realloc(M, nclus * nclus * sizeof(int));

    for (size_t id = 1; id < (nclus + 1); id++)
    {
        for (size_t i = 0; i < n; i++)
        {
            if (id != c[i]) // c[i]: cluster of row i of colCompressed/row
                continue;
            for (size_t j = 0; j < nclus; j++)
            {
                M[(id - 1) * nclus + j] += colCompressed[i][j]; // compress rows of colCompressed by summing the rows of each cluster to the row of M the cluster id points to
            }
        }
    }
};