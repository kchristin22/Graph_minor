#include <iostream>
#include <sys/time.h>
#include <string.h>
#include <fstream>
#include <math.h>
#include <iomanip>
#include <algorithm>
#include "writeMM.hpp"
#include "readMM.hpp"
#include "GMsequential.hpp"
#include "GMopenMP.hpp"
#include "GMpthreads.hpp"

#define NUM_THREADS 4

bool areCSRVectorsEqual(CSR &csrSeq, CSR &csrPar)
{
    if (csrSeq.row != csrPar.row || csrSeq.col.size() != csrPar.col.size() || csrSeq.val.size() != csrPar.val.size() || csrSeq.col.size() != csrSeq.val.size())
        return false;

    for (size_t i = 0; i < csrSeq.row.size(); i++)
    {
        std::vector<size_t>::iterator startRange = csrSeq.col.begin() + csrSeq.row[i];

        for (size_t j = csrSeq.row[i]; j < csrSeq.row[i + 1]; j++)
        {

            std::vector<size_t>::iterator colIndex = std::find(startRange, startRange + csrSeq.row[i + 1], csrPar.col[j]); // find this column in the range of columns of the same row in Seq
                                                                                                                           // (each column appears at most once in each row)
            // If col[j] is not found, vectors are not equal
            if (colIndex == csrSeq.col.end())
            {
                printf("colIndex == colSeq.end()\n");
                return false;
            }

            // Get the index of col[j] in colSeq and thus valSeq
            size_t valIndex = std::distance(csrSeq.col.begin(), colIndex);

            // Check if these values are in the same column cluster
            if (csrPar.val[j] != csrSeq.val[valIndex])
            {
                printf("col[i] != colSeq[colIndex]\n");
                printf(" j = %ld, colIndex = %ld\n", j, valIndex);
                return false;
            }
        }
    }
    // Vectors are equal
    return true;
}

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
    int nz;
    bool symmetric = false;

    verifyMMfile(&Nread, &nz, symmetric, filename);

    size_t nzread = nz;

    std::vector<size_t> conf(Nread, 1);
    char *cfilename = (char *)argv[2];
    std::ifstream cfile(cfilename);

    size_t c, index = 0;
    while (cfile >> c)
    {
        conf[index++] = c;
    }

    /* Parsing the exponential data of the roads_usa cluster ids
    // std::ofstream output;
    // output.open("parsing.txt");

    // size_t sc[10];

    // for (size_t i = 0; i < 5; i++)
    // {
    //     cfile >> std::scientific >> sc[i];
    //     output << std::fixed << std::setprecision(0) << sc[i] << std::endl;
    // }
    // output.close();
    */

    printf("start\n");
    std::vector<size_t> I(nzread, 0);
    std::vector<size_t> J(nzread, 0);
    std::vector<uint64_t> V(nzread, 0);

    std::vector<size_t> row(Nread + 1, 0);
    std::vector<size_t> col(nzread, 0);
    std::vector<uint64_t> val(nzread, 0);

    COO coo = {I, J, V};
    CSR csr = {row, col, val};

    readMM(I, J, V, filename, nzread);
    coo_to_csr(csr, nzread, coo, Nread, symmetric);

    std::vector<size_t> rowM(Nread + 1, 0);
    std::vector<size_t> colM(nzread, 0);
    std::vector<uint64_t> valM(nzread, 0);

    CSR csrM = {rowM, colM, valM};

    struct timeval start, end;
    gettimeofday(&start, NULL);
    seq(csrM, csr, conf);
    gettimeofday(&end, NULL);
    printf("seq time: %ld\n", ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)));

    printf("row size: %ld\n", csrM.row.size());
    printf("col size: %ld\n", csrM.col.size());
    printf("val size: %ld\n", csrM.val.size());

    std::vector<size_t> rowMP(Nread + 1, 0);
    std::vector<size_t> colMP(nzread, 0);
    std::vector<uint64_t> valMP(nzread, 0);

    CSR csrMP = {rowMP, colMP, valMP};

    gettimeofday(&start, NULL);
    GMopenMP(csrMP, csr, conf, numThreads);
    gettimeofday(&end, NULL);
    printf("parallel time: %ld\n", ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)));

    printf("Are seq and openmp vectors equal? %d\n", areCSRVectorsEqual(csrM, csrMP));

    std::vector<size_t> rowPt(Nread + 1, 0);
    std::vector<size_t> colPt(nzread, 0);
    std::vector<uint64_t> valPt(nzread, 0);

    CSR csrPt = {rowPt, colPt, valPt};

    gettimeofday(&start, NULL);
    GMpthreads(csrPt, csr, conf, numThreads);
    gettimeofday(&end, NULL);
    printf("pthread time: %ld\n", ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)));

    printf("Are seq and pthreads vectors equal? %d\n", areCSRVectorsEqual(csrM, csrPt));

    return 0;
}
