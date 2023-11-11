#pragma once

#include "vector"
#include "cstddef"

inline void numClusters(size_t &nclus, std::vector<size_t> &c);

/* Dense matrix implementation */
void seq(std::vector<int> &M, std::vector<int> &A, std::vector<size_t> &c);

/* CSR matrix implementation */
void seq(std::vector<size_t> &rowM, std::vector<size_t> &colM, std::vector<int> &valM,
         std::vector<size_t> &row, std::vector<size_t> &col, std::vector<int> &val, std::vector<size_t> &c);

/* Calculate result using dense matrix implementation for intermediate steps */
void seqDenseCSR(std::vector<size_t> &rowM, std::vector<size_t> &colM, std::vector<int> &valM,
                 std::vector<size_t> &row, std::vector<size_t> &col, std::vector<int> &val, std::vector<size_t> &c);
