
#pragma once

#include "omp.h"
#include "stdio.h"
#include "vector"

inline void numClusters(size_t &nclus, std::vector<size_t> &c);

void openMP(std::vector<size_t> &rowM, std::vector<size_t> &colM, std::vector<int> &valM,
            std::vector<size_t> &row, std::vector<size_t> &col, std::vector<int> &val, std::vector<size_t> &c);