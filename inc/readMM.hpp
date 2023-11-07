#include "mmio.h"
#include "cstdlib"
#include "vector"

void verifyMMfile(int *N, int *nz, char *filename);

void readMM(std::vector<int> *A, char *filename, int N, int nz);