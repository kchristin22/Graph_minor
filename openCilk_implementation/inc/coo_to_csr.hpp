#pragma once

#include <cstdlib>
#include <vector>
#include <stdint.h>
#include <stdio.h>

/* COO representation: I, J, V */
struct COO
{
    std::vector<size_t> &I;
    std::vector<size_t> &J;
    std::vector<uint32_t> &V;
};

/* CSR representation: row, col, val */
struct CSR
{
    std::vector<size_t> &row;
    std::vector<size_t> &col;
    std::vector<uint32_t> &val;
};

/* Receives COO format as input (I, J, V) and transforms it to CSR (row, col, val)
 * @params:
 *    csr (output): CSR representation of the input COO
 *    coo (input): COO representation of the input matrix
 *    N (input): total number of rows/columns of the input matrix / dimensions of the input matrix
 *    symmetric (input): true if the input matrix is symmetric, false otherwise
 */
void coo_to_csr(CSR &csr, const COO &coo, const size_t N, const bool symmetric);