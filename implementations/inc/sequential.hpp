#pragma once

#include "vector"
#include "cstddef"
#include "stdint.h"
#include "coo_to_csr.hpp"

inline void numClusters(size_t &nclus, const std::vector<size_t> &c);

/* Dense matrix implementation */
void seq(std::vector<int> &M, const std::vector<int> &A, const std::vector<size_t> &c);

/* CSR matrix implementation */
void seq(CSR &csrM, const CSR &csr, const std::vector<size_t> &c);

/* Calculate result using dense matrix implementation for intermediate steps */
void seqDenseCSR(std::vector<size_t> &rowM, std::vector<size_t> &colM, std::vector<uint32_t> &valM,
                 const std::vector<size_t> &row, const std::vector<size_t> &col, const std::vector<uint32_t> &val, const std::vector<size_t> &c);
