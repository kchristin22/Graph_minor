#pragma once

#include "stdio.h"
#include "cilk/cilk.h"
#include "vector"
#include "stdint.h"
#include "coo_to_csr.hpp"

#define ELEMENTS_PER_CACHE_LINE_INT (64 / sizeof(int))
#define ELEMENTS_PER_CACHE_LINE_SIZE_T (64 / sizeof(size_t))

inline void numClusters(size_t &nclus, std::vector<size_t> &c);

void GMopenCilk(CSR &csrM, CSR &csr, std::vector<size_t> &c);
