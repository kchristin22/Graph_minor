
#pragma once

#include "omp.h"
#include "stdio.h"
#include "vector"
#include "stdint.h"
#include "coo_to_csr.hpp"

inline void numClusters(size_t &nclus, const std::vector<size_t> &c);

void openMP(CSR &csrM, const CSR &csr, const std::vector<size_t> &c);