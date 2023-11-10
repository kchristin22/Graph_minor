#pragma once

#include "vector"
#include "cstddef"

void seqDense(std::vector<int> &M, std::vector<int> &A, std::vector<size_t> &c); // mixture (to be deleted after testing)

/* Dense matrix implementation */
void seq(std::vector<int> &M, std::vector<int> &A, std::vector<size_t> &c);

/* CSR matrix implementation */
void seq(std::vector<size_t> &rowM, std::vector<size_t> &colM, std::vector<int> &valM,
         std::vector<size_t> &row, std::vector<size_t> &col, std::vector<int> &val, std::vector<size_t> &c);
