#pragma once

#include <cstdlib>
#include <vector>
#include <stdint.h>
#include <stdio.h>

struct COO
{
    std::vector<size_t> &I;
    std::vector<size_t> &J;
    std::vector<uint32_t> &V;
};

struct CSR
{
    std::vector<size_t> &row;
    std::vector<size_t> &col;
    std::vector<uint32_t> &val;
};

/* Receives COO format as input (I, J, V) and transforms it to CSR (row, col, val) */
void coo_to_csr(CSR &csr, const COO &coo, const size_t N, const bool symmetrical, const uint32_t numThreads);