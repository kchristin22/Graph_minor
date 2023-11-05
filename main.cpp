#include "iostream"
#include "stdio.h"
#include "cstdlib"
#include "sequential.hpp"
#include "mmio.h"

int main()
{
    std::vector<std::vector<size_t>> v(1);
    v[0].push_back(1);
    std::vector<size_t> c(1);

    MM_typecode matcode;
    printf("test %d\n", mm_is_valid(matcode));
    seq(v, c, 1);
    return 0;
}
