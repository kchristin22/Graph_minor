#include "mmio.h"
#include "cstdlib"
#include "vector"

void verifyMMfile(int *N, int *nz, char *filename);

void readMM(std::vector<int> &A, char *filename, int N, int nz);

void readMM(std::vector<size_t> &I, std::vector<size_t> &J, std::vector<int> &V, char *filename, int N, int nz);