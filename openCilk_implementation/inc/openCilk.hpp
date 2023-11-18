#pragma once

#include "stdio.h"
#include "cilk/cilk.h"
#include "vector"
#include "stdint.h"

#define ELEMENTS_PER_CACHE_LINE_INT (64 / sizeof(int))
#define ELEMENTS_PER_CACHE_LINE_SIZE_T (64 / sizeof(size_t))

inline void numClusters(size_t &nclus, std::vector<size_t> &c);

void GMopenCilk(std::vector<size_t> &rowM, std::vector<size_t> &colM, std::vector<uint32_t> &valM,
                std::vector<size_t> &row, std::vector<size_t> &col, std::vector<uint32_t> &val, std::vector<size_t> &c);
