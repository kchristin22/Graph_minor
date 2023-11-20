#include <atomic>
#include <omp.h>
#include "GMopenMP.hpp"


void coo_to_csr(CSR &csr, const COO &coo, const size_t N, const bool symmetrical, const uint32_t numThreads)
{

    if (coo.I.size() != coo.J.size() || coo.I.size() != coo.V.size())
    {
        printf("Error: at least one of the pairs has unmatched dimensions: (I,J), (I,V) \n");
        exit(1);
    }
    if ((csr.row.size() != (N + 1)) || (csr.col.size() != coo.J.size()) || (csr.val.size() != coo.V.size()))
    {
        printf("Error: at least one of the outputs has wrong dimensions \n");
        exit(1);
    }

    size_t nz = csr.val.size();

    size_t x = 0, rowIndex = 0, localCount;
    std::atomic<size_t> count(0);

    if (symmetrical) // we can take advantage of the J vector being sorted and use that vector as the row index
    {
        csr.row[0] = 0;
        for (size_t index = 0; index < nz; index++)
        {
            if (coo.J[index] != x)
            {
                csr.row[++rowIndex] = index; // offset of the next row
                x = coo.J[index];            // new row
            }
            csr.col[index] = coo.I[index];
            csr.val[index] = coo.V[index];
        }
    }
    else
    {
        size_t chunk;
        calChunk(chunk, N, 16, numThreads);

        for (size_t index = 0; index < N; index++)
        {
            csr.row[index] = count;

#pragma omp parallel num_threads(numThreads)
#pragma omp for nowait private(localCount) schedule(dynamic, chunk)
            for (size_t j = 0; j < nz; j++)
            {
                if (coo.I[j] != index)
                    continue;

                localCount = count.fetch_add(1);
                csr.col[localCount] = coo.J[j];
                csr.val[localCount] = coo.V[j];
            }
        }

        csr.row[N] = count;
    }
}