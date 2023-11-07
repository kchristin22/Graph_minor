#include "writeMM.hpp"
#include "cstdlib"
#include "algorithm"

void writeMM(char *filename, int N, int nz)
{
  MM_typecode matcode;

  mm_initialize_typecode(&matcode);
  mm_set_matrix(&matcode);
  mm_set_coordinate(&matcode);
  mm_set_real(&matcode);

  int i, I[nz], J[nz], val[nz];
  for (i = 0; i < nz; i++)
  {
    I[i] = rand() % N;
    J[i] = rand() % N;
    val[i] = rand() % 10;
  }

  std::sort(I, I + nz);

  FILE *f = fopen(filename, "w");

  mm_write_banner(f, matcode);
  mm_write_mtx_crd_size(f, N, N, nz);

  /* NOTE: matrix market files use 1-based indices, i.e. first element
of a vector has index 1, not 0. */

  for (i = 0; i < nz; i++)
    fprintf(f, "%d %d %d\n", I[i] + 1, J[i] + 1, val[i]);

  fclose(f);
}
