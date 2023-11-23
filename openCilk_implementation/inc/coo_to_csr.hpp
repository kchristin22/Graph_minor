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

/* Receives COO format as input (I, J, V) and transforms it to CSR (row, col, val) (from scipy library)
 * @params:
 *    n_row (input): total number of rows of the input matrix
 *    n_col (input): total number of columns of the input matrix
 *    nnz (input): total number of non-zero elements of the input matrix
 *    Ai (input): row indexes of the non-zero elements of the input matrix
 *    Aj (input): column indexes of the non-zero elements of the input matrix
 *    Ax (input): values of the non-zero elements of the input matrix
 *    Bp (output): CSR row vector of the input matrix
 *    Bj (output): CSR col vector of the input matrix
 *    Bx (output): CSR val vector of the input matrix
 */
template <class I, class T>
void coo_tocsr(const I n_row,
               const I n_col,
               const I nnz,
               const I Ai[],
               const I Aj[],
               const T Ax[],
               I Bp[],
               I Bj[],
               T Bx[])
{
    // compute number of non-zero entries per row of A
    std::fill(Bp, Bp + n_row, 0);

    for (I n = 0; n < nnz; n++)
    {
        Bp[Ai[n]]++;
    }

    // cumsum the nnz per row to get Bp[]
    for (I i = 0, cumsum = 0; i < n_row; i++)
    {
        I temp = Bp[i];
        Bp[i] = cumsum;
        cumsum += temp;
    }
    Bp[n_row] = nnz;

    // write Aj,Ax into Bj,Bx
    for (I n = 0; n < nnz; n++)
    {
        I row = Ai[n];
        I dest = Bp[row];

        Bj[dest] = Aj[n];
        Bx[dest] = Ax[n];

        Bp[row]++;
    }

    for (I i = 0, last = 0; i <= n_row; i++)
    {
        I temp = Bp[i];
        Bp[i] = last;
        last = temp;
    }

    // now Bp,Bj,Bx form a CSR representation (with possible duplicates)
}

/* Receives COO format as input (I, J, V) and transforms it to CSR (row, col, val)
 * @params:
 *    csr (output): CSR representation of the input COO
 *    coo (input): COO representation of the input matrix
 *    N (input): total number of rows/columns of the input matrix / dimensions of the input matrix
 */
void coo_to_csr(CSR &csr, const COO &coo, const size_t n);

/* Receives COO format as input (I, J, V) and transforms it to CSR (row, col, val)
 * @params:
 *    csr (output): CSR representation of the input COO
 *    coo (input): COO representation of the input matrix
 *    N (input): total number of rows/columns of the input matrix / dimensions of the input matrix
 *    symmetric (input): true if the input matrix is symmetric, false otherwise
 */
void slow_coo_to_csr(CSR &csr, const COO &coo, const size_t N, const bool symmetric);