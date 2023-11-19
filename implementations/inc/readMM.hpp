#include "mmio.h"
#include "stdlib.h"
#include "vector"
#include "stdint.h"

void verifyMMfile(int *N, int *nz, const char *filename);

void readMM(std::vector<int> &A, const char *filename, const int N, const int nz);

void readMM(std::vector<size_t> &I, std::vector<size_t> &J, std::vector<uint32_t> &V, const char *filename, const int N, const int nz);