#include "iostream"
#include "stdio.h"
#include "cstdlib"
#include "sequential.hpp"
#include "ctype.h"
#include "cmath"

void seq(std::vector<std::vector<size_t>> A, std::vector<size_t> c, size_t n)
{
    printf("Hello from seq\n");
    // printf("A[0]: %ld\n", A[0].size());
    // printf("c: %ld\n", c.size());
    // printf("size vector: %ld\n", A.size());
    // printf("size ^ 2: %ld\n", (size_t)pow(size, 2));

    if ((A[0].size() != n) || (A.size() != (size_t)(pow(n, 2))) || (c.size() != n))
    {
        printf("Error: Incorrect size of A and c\n");
        exit(1);
    }

    std::vector<size_t> row(A[0].size()); // offset of val vector for each row
    std::vector<size_t> col(A.size());    // col index of each element of val vector
    std::vector<size_t> val(A.size());    // vector containing the non-zero elements of A

    size_t count = 0, index, cond, x;

    for (size_t i = 0; i < n; i++)
    {

        for (int j = 0; j < n; j++)
        {
            x = A[i][j];
            if (x == 0)
                continue;
            count++;
            val[count] = x;
            col[count] = j;
            row[count] = i;
        }
    }
}
