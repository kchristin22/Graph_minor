#pragma once

#include "vector"
#include "cstddef"

void seqDense(std::vector<int> &M, size_t *L, std::vector<int> &A, std::vector<size_t> &c); // mixture (to be deleted after testing)

/* Dense matrix implementation */
void seq(std::vector<int> &M, size_t *L, std::vector<int> &A, std::vector<size_t> &c);
