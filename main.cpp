#include "sequential.hpp"
#include "writeMM.hpp"
#include "readMM.hpp"
#include "openMP.hpp"
#include "cmath"
#include "sys/time.h"

/* Receives COO format as input (I, J, V) and transforms it to CSR (row, col, val) */
void csr(std::vector<size_t> &row, std::vector<size_t> &col, std::vector<int> &val,
         std::vector<size_t> &I, std::vector<size_t> &J, std::vector<int> &V, size_t N, bool symmetrical)
{
    if (I.size() != J.size() || I.size() != V.size())
    {
        printf("Error: at least one of the pairs has unmatched dimensions: (I,J), (I,V) \n");
        exit(1);
    }
    if ((row.size() != N) || (col.size() != J.size()) || (val.size() != V.size()))
    {
        printf("Error: at least one of the outputs has wrong dimensions \n");
        exit(1);
    }

    size_t nz = val.size();

    size_t x = 0, rowIndex = 0, count = 0;

    if (symmetrical) // we can take advantage of the J vector being sorted and use that vector as the row index
    {
        row[0] = 0;
        for (size_t index = 0; index < nz; index++)
        {
            if (J[index] != x)
            {
                row[++rowIndex] = index;
                x = J[index]; // new row
            }
            col[index] = I[index];
            val[index] = V[index];
        }
    }
    else
    {
        size_t chunk = N / 16;
        size_t numThreads = 4;
        if (!chunk)
        {
            chunk = N;
            numThreads = 1;
        }

        // #pragma omp parallel for reduction(+ : count) schedule(dynamic, chunk) num_threads(numThreads)
        for (size_t index = 0; index < N; index++)
        {
            row[index] = count;
            for (size_t j = 0; j < nz; j++)
            {
                if (I[j] != index)
                    continue;
                col[count] = J[j];
                val[count] = V[j];
                count++;
            }
        }
    }
}

int main(int argc, char *argv[])
{

    if (argc < 2)
    {
        fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
        exit(1);
    }

    char *filename = (char *)argv[1];

    int Nread;
    int nzread;

    verifyMMfile(&Nread, &nzread, filename);

    // std::vector<int> A(Nread * Nread, 0);

    // readMM(A, filename, Nread, nzread);

    std::vector<size_t> conf(Nread);
    for (int i = 0; i < Nread; i++)
    {
        conf[i] = rand() % 24 + 1;
        //     printf("%ld ", conf[i]);
    }
    // printf("\n");

    // std::vector<int> M(Nread * Nread, 0); // nz x nz max

    // struct timeval before, after;
    // gettimeofday(&before, NULL);
    // seq(M, A, conf);
    // gettimeofday(&after, NULL);
    // printf("dense time: %ld\n", ((after.tv_sec * 1000000 + after.tv_usec) - (before.tv_sec * 1000000 + before.tv_usec)));
    // A.~vector();

    // size_t L = (size_t)sqrt(M.size()); // maybe dangerous, consider adding it again as an input arg

    // for (size_t i = 0; i < L; i++)
    // {
    //     for (size_t j = 0; j < L; j++)
    //     {
    //         printf("%d ", M[(i * L) + j]);
    //     }
    //     printf("\n");
    // }
    std::vector<size_t> I(nzread, 0);
    std::vector<size_t> J(nzread, 0);
    std::vector<int> V(nzread, 0);

    std::vector<size_t> row(Nread, 0);
    std::vector<size_t> col(nzread, 0);
    std::vector<int> val(nzread, 0);

    readMM(I, J, V, filename, Nread, nzread);
    csr(row, col, val, I, J, V, Nread, false);

    // for (size_t i = 0; i < nzread; i++)
    // {
    //     printf("%d %d %d \n", I[i], J[i], V[i]);
    // }
    // printf("CSR:\n");
    // printf("row: ");
    // for (size_t i = 0; i < row.size(); i++)
    // {
    //     printf("%d ", row[i]);
    // }
    // printf("\ncol: ");
    // for (size_t i = 0; i < col.size(); i++)
    // {
    //     printf("%d ", col[i]);
    // }
    // printf("\nval: ");
    // for (size_t i = 0; i < val.size(); i++)
    // {
    //     printf("%d ", val[i]);
    // }
    // printf("\n");

    std::vector<size_t> rowM(Nread, 0);
    std::vector<size_t> colM(nzread, 0);
    std::vector<int> valM(nzread, 0);

    struct timeval start, end;
    gettimeofday(&start, NULL);
    seq(rowM, colM, valM, row, col, val, conf);
    gettimeofday(&end, NULL);
    printf("seq time: %ld\n", ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)));

    // printf("rowM: ");
    // for (size_t i = 0; i < rowM.size(); i++)
    // {
    //     printf("%ld ", rowM[i]);
    // }
    // printf("\ncolM: ");
    // for (size_t i = 0; i < colM.size(); i++)
    // {
    //     printf("%ld ", colM[i]);
    // }
    // printf("\nvalM: ");
    // for (size_t i = 0; i < valM.size(); i++)
    // {
    //     printf("%d ", valM[i]);
    // }
    // printf("\n");

    colM.resize(nzread, 0);
    valM.resize(nzread, 0);
    rowM.resize(Nread, 0);

    gettimeofday(&start, NULL);
    openMP(rowM, colM, valM, row, col, val, conf);
    gettimeofday(&end, NULL);
    printf("parallel time: %ld\n", ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)));

    // printf("rowM: ");
    // for (size_t i = 0; i < rowM.size(); i++)
    // {
    //     printf("%ld ", rowM[i]);
    // }
    // printf("\ncolM: ");
    // for (size_t i = 0; i < colM.size(); i++)
    // {
    //     printf("%ld ", colM[i]);
    // }
    // printf("\nvalM: ");
    // for (size_t i = 0; i < valM.size(); i++)
    // {
    //     printf("%d ", valM[i]);
    // }
    // printf("\n");

    return 0;
}
