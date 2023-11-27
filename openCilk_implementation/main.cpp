#include <sys/time.h>
#include <fstream>
#include "writeMM.hpp"
#include "readMM.hpp"
#include "GMopenCilk.hpp"

int main(int argc, char *argv[])
{

    if (argc < 2)
    {
        fprintf(stderr, "Usage: %s ${martix-market-filename} ${c-matrix-filename}\n", argv[0]);
        exit(1);
    }

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

    std::vector<size_t> I(nzread, 0);
    std::vector<size_t> J(nzread, 0);
    std::vector<uint64_t> V(nzread, 0);

    std::vector<size_t> row(Nread + 1, 0);
    std::vector<size_t> col(nzread, 0);
    std::vector<uint64_t> val(nzread, 0);

    CSR csr = {row, col, val};
    COO coo = {I, J, V};

    readMM(I, J, V, filename, nzread);

    coo_to_csr(csr, nzread, coo, Nread, symmetric);

    std::vector<size_t> rowM(Nread + 1, 0);
    std::vector<size_t> colM(nzread, 0);
    std::vector<uint64_t> valM(nzread, 0);

    CSR csrM = {rowM, colM, valM};

    struct timeval start, end;

    gettimeofday(&start, NULL);
    GMopenCilk(csrM, csr, conf);
    gettimeofday(&end, NULL);
    printf("openCilk time: %ld\n", ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)));

    printf("row size: %ld\n", csrM.row.size());
    printf("col size: %ld\n", csrM.col.size());
    printf("val size: %ld\n", csrM.val.size());


    return 0;
}
