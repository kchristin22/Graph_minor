#include <iostream>
#include <sys/time.h>
#include <string.h>
#include <fstream>
#include "writeMM.hpp"
#include "readMM.hpp"
#include "GMsequential.hpp"
#include "GMopenMP.hpp"
#include "GMpthreads.hpp"

#define NUM_THREADS 4

int main(int argc, char *argv[])
{

    uint32_t numThreads = NUM_THREADS;

    if (argc < 3)
    {
        fprintf(stderr, "Usage: %s ${martix-market-filename} ${c-matrix-filename} ${number-of-threads}\n", argv[0]);
        exit(1);
    }
    else if (argc == 4)
        if (atoi(argv[3]) > 0)
            numThreads = atoi(argv[3]);

    std::cout << "Running with numThreads: " << numThreads << std::endl;

    char *filename = (char *)argv[1];

    int Nread;
    int nzread;

    verifyMMfile(&Nread, &nzread, filename);

    std::vector<size_t> conf(Nread, 1);
    char *cfilename = (char *)argv[2];
    std::ifstream cfile(cfilename);

    size_t c, index = 0;
    while (cfile >> c)
    {
        conf[index++] = c;
    }

    // printf("conf: ");
    // for (size_t i = 0; i < conf.size(); i++)
    // {
    //     printf("%ld ", conf[i]);
    // }

    // for (int i = 0; i < 25000; i++)
    // {
    //     conf[i] = i + 1;
    //     // printf("%ld ", conf[i]);
    // }
    // printf("\n");

    printf("start\n");
    std::vector<size_t> I(nzread, 0);
    std::vector<size_t> J(nzread, 0);
    std::vector<uint32_t> V(nzread, 0);

    std::vector<size_t> row(Nread + 1, 0);
    std::vector<size_t> col(nzread, 0);
    std::vector<uint32_t> val(nzread, 0);

    COO coo = {I, J, V};
    CSR csr = {row, col, val};

    readMM(I, J, V, filename, nzread);
    coo_to_csr(csr, coo, Nread, false, numThreads);

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

    printf("row size: %ld\n", csrM.row.size());
    printf("col size: %ld\n", csrM.col.size());
    printf("val size: %ld\n", csrM.val.size());

    std::vector<size_t> rowM2(Nread + 1, 0);
    rowM2 = rowM;
    std::vector<size_t> colM2(nzread, 0);
    colM2 = colM;

    colM.resize(nzread, 0);
    valM.resize(nzread, 0);
    rowM.resize(Nread + 1, 0);

    gettimeofday(&start, NULL);
    GMopenMP(csrM, csr, conf, numThreads);
    gettimeofday(&end, NULL);
    printf("parallel time: %ld\n", ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)));

    std::vector<size_t> rowM3(Nread + 1, 0);
    rowM3 = rowM;
    std::vector<size_t> colM3(nzread, 0);
    colM3 = colM;

    if (rowM2 != rowM3 || colM2.size() != colM3.size())
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

    gettimeofday(&start, NULL); // move this inside?
    GMpthreads(csrM, csr, conf, numThreads);
    gettimeofday(&end, NULL);
    printf("pthread time: %ld\n", ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)));
    std::vector<size_t> rowM4(Nread + 1, 0);
    rowM4 = rowM;
    std::vector<size_t> colM4(nzread, 0);
    colM4 = colM;

    if (rowM2 != rowM4 || colM2.size() != colM4.size()){
        printf("seq not in accordance to pthreads\n");
        printf("rowM2.size: %ld, rowM4.size: %ld\n", rowM2.size(), rowM4.size());
        printf("colM2.size: %ld, colM4.size: %ld\n", colM2.size(), colM4.size());
    }



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

    return 0;
}
