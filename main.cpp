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

    std::vector<int> v(16);
    v = {0, 2, 3, 4, 5, 0, 2, 2, 0, 0, 0, 1, 0, 9, 0, 8};
    std::vector<size_t> c(4);
    c = {1, 2, 1, 2};

    std::vector<int> M(16, 0);
    size_t L;

    seq(&M, &L, &v, &c);
    for (size_t i = 0; i < L; i++)
    {
        for (size_t j = 0; j < L; j++)
        {
            printf("%d ", M[(i * L) + j]);
        }
        printf("\n");
    }

    char *filename = (char *)argv[1];

    int Nread;
    int nzread;

    writeMM(filename, 20, 10);
    verifyMMfile(&Nread, &nzread, filename);
    std::vector<int> A(Nread * Nread, 0);
    readMM(&A, filename, Nread, nzread);

    printMP();
    return 0;
}
