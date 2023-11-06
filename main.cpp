#include "iostream"
#include "stdio.h"
#include "cstdlib"
#include "sequential.hpp"
#include "mmio.h"

int main()
{
    size_t v[4][4] = {{0, 2, 3, 4}, {5, 0, 2, 2}, {0, 0, 0, 1}, {0, 9, 0, 8}};
    size_t c[4] = {1, 2, 1, 2};
    MM_typecode matcode;
    printf("test %d\n", mm_is_valid(matcode));
    int *M = new int[16]{0};
    size_t L;
    seq<4>(M, &L, v, c);
    for (size_t i = 0; i < L; i++)
    {
        for (size_t j = 0; j < L; j++)
        {
            printf("%d ", M[i * L + j]);
        }
        printf("\n");
    }
    return 0;
}
