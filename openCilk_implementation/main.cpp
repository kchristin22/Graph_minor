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

    // for (int i = 0; i < 8000; i++)
    // {
    //     conf[i] = i + 1;
    //     // printf("%ld ", conf[i]);
    // }
    // printf("\n");

    std::vector<size_t> I(nzread, 0);
    std::vector<size_t> J(nzread, 0);
    std::vector<uint32_t> V(nzread, 0);

    std::vector<size_t> row(Nread + 1, 0);
    std::vector<size_t> col(nzread, 0);
    std::vector<uint32_t> val(nzread, 0);

    CSR csr = {row, col, val};
    COO coo = {I, J, V};

    readMM(I, J, V, filename, nzread);

    coo_to_csr(csr, coo, Nread, false);

    // for (size_t i = 0; i < nzread; i++)
    // {
    //     printf("%d %d %d \n", I[i], J[i], V[i]);
    // }

    // printf("CSR:\n");
    // printf("row: ");
    // for (size_t i = 0; i < csr.row.size(); i++)
    // {
    //     printf("%ld ", csr.row[i]);
    // }
    // printf("\ncol: ");
    // for (size_t i = 0; i < col.size(); i++)
    // {
    //     printf("%ld ", csr.col[i]);
    // }
    // printf("\nval: ");
    // for (size_t i = 0; i < val.size(); i++)
    // {
    //     printf("%d ", csr.val[i]);
    // }
    // printf("\n");

    std::vector<size_t> rowM(Nread + 1, 0);
    std::vector<size_t> colM(nzread, 0);
    std::vector<uint32_t> valM(nzread, 0);

    CSR csrM = {rowM, colM, valM};

    struct timeval start, end;

    gettimeofday(&start, NULL);
    GMopenCilk(csrM, csr, conf);
    gettimeofday(&end, NULL);
    printf("openCilk time: %ld\n", ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)));

    printf("row size: %ld\n", csrM.row.size());
    printf("col size: %ld\n", csrM.col.size());
    printf("val size: %ld\n", csrM.val.size());
    // printf("rowM: ");
    // for (size_t i = 0; i < csrM.row.size(); i++)
    // {
    //     printf("%ld ", csrM.row[i]);
    // }
    // printf("\ncolM: ");
    // for (size_t i = 0; i < csrM.col.size(); i++)
    // {
    //     printf("%ld ", csrM.col[i]);
    // }
    // printf("\nvalM: ");
    // for (size_t i = 0; i < csrM.val.size(); i++)
    // {
    //     printf("%d ", csrM.val[i]);
    // }
    // printf("\n");

    return 0;
}
