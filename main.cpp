#include "sequential.hpp"
#include "writeMM.hpp"
#include "readMM.hpp"
#include "openMP.hpp"

#define N 20
#define nz 10

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

    writeMM(filename, 20, 10);
    verifyMMfile(&Nread, &nzread, filename);

    std::vector<int> A(Nread * Nread, 0);
    readMM(A, filename, Nread, nzread);

    std::vector<size_t> conf(Nread);
    for (size_t i = 0; i < Nread; i++)
    {
        conf[i] = rand() % 5;
        printf("%ld ", conf[i]);
    }
    printf("\n");

    std::vector<int> M(Nread * Nread, 0); // nz x nz max
    size_t L;

    seq(M, &L, A, conf);
    for (size_t i = 0; i < L; i++)
    {
        for (size_t j = 0; j < L; j++)
        {
            printf("%d ", M[(i * L) + j]);
        }
        printf("\n");
    }

    printMP();
    return 0;
}
