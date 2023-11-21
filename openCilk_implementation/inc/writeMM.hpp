#include "mmio.h"

/* Write a matrixof certain dimensions and number of non-zero elements in a MM file with
 * @params:
 *   filename (input): name of the file to be written
 *   N (input): dimension of the matrix to be written
 *   nz (input): number of non-zero elements of the matrix to be written
 */
void writeMM(const char *filename, const int N, const int nz);