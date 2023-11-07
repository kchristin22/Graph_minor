#include "iostream"
#include "stdio.h"
#include "cstdlib"
#include "sequential.hpp"
#include "writeMM.hpp"
// #include "readMM.hpp
#include "openMP.hpp"

#define nz 10
#define N 20

int main()
{
    size_t v[4][4] = {{0, 2, 3, 4}, {5, 0, 2, 2}, {0, 0, 0, 1}, {0, 9, 0, 8}};
    size_t c[4] = {1, 2, 1, 2};

    std::vector<int> M(16, 0);
    size_t L;

    seq<4>(&M, &L, v, c);
    for (size_t i = 0; i < L; i++)
    {
        for (size_t j = 0; j < L; j++)
        {
            printf("%d ", M[(i * L) + j]);
        }
        printf("\n");
    }

    // writeMM(nz, N);

    // check in main:
    // if ((A[0].size() != n) || (A.size() != (size_t)(pow(n, 2))) || (c.size() != n))
    // {
    //     printf("Error: Incorrect size of A or c\n");
    //     exit(1);
    // }

    printMP();
    return 0;
}
