#include <atomic>
#include <omp.h>
#include <iostream>
#include "GMopenMP.hpp"

void coo_to_csr(CSR &csr, size_t &nnz, const COO &coo, const size_t n, const bool symmetric)
{

    nnz = coo.I.size();

    if (coo.I.size() != coo.J.size() || coo.I.size() != coo.V.size())
    {
        std::cout << "coo_to_csr: I, J, V must have the same size" << std::endl;
        exit(1);
    }

    if (symmetric)
    {
        size_t allCount = 0;
        for (size_t i = 0; i < nnz; i++)
        {
            if (coo.I[i] != coo.J[i])
            {
                coo.I.push_back(coo.J[i]);
                coo.J.push_back(coo.I[i]);
                coo.V.push_back(coo.V[i]);
                allCount++;
            }
        }
        nnz += allCount;
    }

    csr.row.resize(n + 1);
    csr.col.resize(nnz);
    csr.val.resize(nnz);

    coo_tocsr(n, n, nnz, coo.I.data(), coo.J.data(), coo.V.data(), csr.row.data(), csr.col.data(), csr.val.data());
}

void slow_coo_to_csr(CSR &csr, const COO &coo, const size_t N, const bool symmetric, const uint32_t numThreads)
{
    // the input vectors must have the same size by definition (COO format)
    if (coo.I.size() != coo.J.size() || coo.I.size() != coo.V.size())
    {
        printf("Error: at least one of the pairs has unmatched dimensions: (I,J), (I,V) \n");
        exit(1);
    }
    // the corresponding output vectors of J and V must have the same size as them,
    // and the row vector should be of size equal to the rows of the matrix, plus one (the latter is an implementation specific)
    if ((csr.row.size() != (N + 1)) || (csr.col.size() != coo.J.size()) || (csr.val.size() != coo.V.size()))
    {
        printf("Error: at least one of the outputs has wrong dimensions \n");
        exit(1);
    }

    size_t nz = csr.val.size(); // number of non-zero elements

    size_t x = 0, rowIndex = 0;

    if (symmetric)      // we can take advantage of the J vector being sorted and use that vector as the row index (swap I and J)
    {                   // this part has not been tested
        csr.row[0] = 0; // or csr.row[0] = x;
        for (size_t index = 0; index < nz; index++)
        {
            if (coo.J[index] != x) // the row index has changed
            {
                csr.row[++rowIndex] = index; // offset of the col and val vectors of the next row
                x = coo.J[index];            // update the row index
            }
            csr.col[index] = coo.I[index];
            csr.val[index] = coo.V[index];
        }
    }
    else
    {
        std::atomic<size_t> count(0); // keep track of the number of non-zero elements in each row
                                      // this variable is atomic because it is shared between threads and it is updated by each thread
        size_t localCount;            // private variable of each thread used for storing the value of count fetched at this point before updating it

        size_t chunk;
        calChunk(chunk, N, ELEMENTS_PER_CACHE_LINE_INT, numThreads); // calculate the chunk size for the parallel region to minimize false sharing

        for (size_t index = 0; index < N; index++)
        {
            csr.row[index] = count;

#pragma omp parallel num_threads(numThreads)
#pragma omp for nowait private(localCount) schedule(dynamic, chunk)
            for (size_t j = 0; j < nz; j++)
            {
                if (coo.I[j] != index)
                    continue;

                localCount = count.fetch_add(1); // the variable is updated only upon reading and storing the previous value,
                                                 // to avoid skipping index positions
                csr.col[localCount] = coo.J[j];
                csr.val[localCount] = coo.V[j];
            }
        }

        csr.row[N] = count; // this element is only used to store the range of the col and val vectors corresponding to the last row
    }
}