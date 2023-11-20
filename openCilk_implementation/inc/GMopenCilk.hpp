#pragma once

#include <cilk/cilk.h>
#include "coo_to_csr.hpp"

#define ELEMENTS_PER_CACHE_LINE_INT (64 / sizeof(int))
#define ELEMENTS_PER_CACHE_LINE_SIZE_T (64 / sizeof(size_t))

inline void numClusters(size_t &nclus, const std::vector<size_t> &c, const uint32_t numThreads);

void GMopenCilk(CSR &csrM, const CSR &csr, const std::vector<size_t> &c, const uint32_t numThreads);
