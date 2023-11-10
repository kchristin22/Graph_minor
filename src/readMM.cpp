#include "readMM.hpp"

void verifyMMfile(int *N, int *nz, char *filename)
{
    FILE *f;
    if ((f = fopen(filename, "r")) == NULL)
    {
        printf("Error: Cannot open file %s\n", filename);
        exit(1);
    }

    MM_typecode matcode;
    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && mm_is_sparse(matcode))
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    if (mm_read_mtx_crd_size(f, N, N, nz) != 0)
    {
        printf("Error: Cannot read matrix size.\n");
        exit(1);
    }

    fclose(f);
}

void readMM(std::vector<int> &A, char *filename, int N, int nz)
{
    FILE *f = fopen(filename, "r");

    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    for (int i = 0; i < 2; i++)
    {
        fscanf(f, "%*[^\n]\n"); // skip first 2 lines, because fscanf doesn't move to the next line if not successfull
    }

    int index, i, j, val;
    for (index = 0; index < (nz); index++)
    {
        fscanf(f, "%d %d %d\n", &i, &j, &val);
        A[(--i) * N + (--j)] = val; // we assume that all items of A are initialized to zero and that the matrix indexes in the file start from 1
    }

    if (f != stdin)
        fclose(f);

    // /************************/
    // /* now write out matrix */
    // /************************/

    // mm_write_mtx_crd_size(stdout, N, N, nz);
    // for (i = 0; i < N; i++)
    // {
    //     for (j = 0; j < N; j++)
    //         fprintf(stdout, "%d ", (*A)[i * N + j]);
    //     fprintf(stdout, "\n");
    // }
}

void readMM(std::vector<size_t> &I, std::vector<size_t> &J, std::vector<int> &V, char *filename, int N, int nz)
{

    FILE *f = fopen(filename, "r");

    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    for (int i = 0; i < 2; i++)
    {
        fscanf(f, "%*[^\n]\n"); // skip first 2 lines, because fscanf doesn't move to the next line if not successfull
    }

    int index;
    for (index = 0; index < (nz); index++)
    {
        fscanf(f, "%d %d %d\n", &I[index], &J[index], &V[index]);
        I[index]--;
        J[index]--;
    }

    if (f != stdin)
        fclose(f);
}
