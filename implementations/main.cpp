#include "GMsequential.hpp"
#include "writeMM.hpp"
#include "readMM.hpp"
#include "GMopenMP.hpp"
#include "GMpthreads.hpp"
#include <sys/time.h>
#include <string.h>

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

    std::vector<size_t> conf(Nread, 1);
    for (int i = 0; i < 1000; i++)
    {
        conf[i] = i + 1;
        // printf("%ld ", conf[i]);
    }
    // printf("\n");

    std::vector<size_t> I(nzread, 0);
    std::vector<size_t> J(nzread, 0);
    std::vector<uint32_t> V(nzread, 0);

    std::vector<size_t> row(Nread + 1, 0);
    std::vector<size_t> col(nzread, 0);
    std::vector<uint32_t> val(nzread, 0);

    COO coo = {I, J, V};
    CSR csr = {row, col, val};

    readMM(I, J, V, filename, Nread, nzread);
    coo_to_csr(csr, coo, Nread, false);

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

    std::vector<size_t> rowM(Nread + 1, 0);
    std::vector<size_t> colM(nzread, 0);
    std::vector<uint32_t> valM(nzread, 0);

    CSR csrM = {rowM, colM, valM};

    struct timeval start, end;
    gettimeofday(&start, NULL);
    seq(csrM, csr, conf);
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
    //     printf("%ld ", valM[i]);
    // }
    // printf("\n");

    std::vector<size_t> rowM2(Nread + 1, 0);
    rowM2 = rowM;

    colM.resize(nzread, 0);
    valM.resize(nzread, 0);
    rowM.resize(Nread + 1, 0);

    gettimeofday(&start, NULL);
    GMopenMP(csrM, csr, conf);
    gettimeofday(&end, NULL);
    printf("parallel time: %ld\n", ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)));

    std::vector<size_t> rowM3(Nread + 1, 0);
    rowM3 = rowM;

    if (rowM2 != rowM3)
        printf("seq not in accordance to openmp\n");
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
    //     printf("%ld ", valM[i]);
    // }
    // printf("\n");

    colM.resize(nzread, 0);
    valM.resize(nzread, 0);
    rowM.resize(Nread + 1, 0);

    gettimeofday(&start, NULL);
    GMpthreads(csrM, csr, conf);
    gettimeofday(&end, NULL);
    printf("pthread time: %ld\n", ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)));
    std::vector<size_t> rowM4(Nread + 1, 0);
    rowM4 = rowM;

    if (rowM2 != rowM4)
        printf("seq not in accordance to pthreads\n");

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
    //     printf("%ld ", valM[i]);
    // }
    printf("\n");

    return 0;
}
