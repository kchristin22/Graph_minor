#pragma once

#include <stdlib.h>
#include <vector>
#include <stdint.h>
#include "mmio.h"

/* Verify the MM file and read the properties of the matrix included
 * @params:
 *    N (output): total number of rows/columns of the input matrix
 *    nz (output): total number of non-zero elements of the input matrix
 *    symmetric (output): true if the input matrix is symmetric, false otherwise
 *    filename (input): name of the file to be read
 */
void verifyMMfile(int *N, int *nz, bool &symmetric, const char *filename);

/* Read the input matrix in COO format from the MM file and store it in an output matrix
 * @params:
 *   A (output): dense matrix representation of the input matrix
 *   filename (input): name of the file to be read
 *   N (input): total number of rows/columns of the input matrix to use for the allocation of the output matrix
 *   nz (input): total number of non-zero elements of the input matrix to use for scanning the input matrix
 */
void readMM(std::vector<int> &A, const char *filename, const int N, const int nz);

/* Read the input matrix in COO format from the MM file and store the COO vectors in output vectors
 * @params:
 *   I (output): row indexes of the non-zero elements of the input matrix
 *   J (output): column indexes of the non-zero elements of the input matrix
 *   V (output): values of the non-zero elements of the input matrix
 *   filename (input): name of the file to be read
 *   nz (input): total number of non-zero elements of the input matrix to use for scanning the input matrix and for the allocation of the output vectors
 */
void readMM(std::vector<size_t> &I, std::vector<size_t> &J, std::vector<uint64_t> &V, const char *filename, const int nz);