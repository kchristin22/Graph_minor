
#pragma once

#include "omp.h"
#include "stdio.h"
#include "vector"
#include "stdint.h"
#include "coo_to_csr.hpp"

inline void numClusters(size_t &nclus, std::vector<size_t> &c);

void openMP(CSR &csrM, CSR &csr, std::vector<size_t> &c);