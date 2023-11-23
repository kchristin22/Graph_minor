#include <atomic>
#include <cilk/opadd_reducer.h>
#include <cilk/cilk_api.h>
#include "GMopenCilk.hpp"

inline void zero_s(void *v) { *(size_t *)v = 0; } // initializer for the reducer
inline void max_s(void *l, void *r)               // operator for the reducer
{
    if (*(size_t *)l < *(size_t *)r)
        *(size_t *)l = *(size_t *)r;
}

inline void numClusters(size_t &nclus, const std::vector<size_t> &c)
{
    size_t n = c.size();
    size_t cilk_reducer(zero_s, max_s) max = 0; // variable to which the max reduction will be applied

    cilk_for(size_t i = 0; i < n; i++)
    {
        max = (c[i] > max) ? c[i] : max;
    }

    nclus = max;
}

// void GMopenCilk(CSR &csrM, const CSR &csr, const std::vector<size_t> &c)
// {
//     size_t n = c.size();
//     size_t nz = csr.val.size();

//     // the row vector should be of size equal to the rows of the matrix, plus one (the latter is an implementation specific)
//     // we assume that the c vector's size is set correctly equal to the number of nodes of the input grapg / rows of the input matrix
//     if (csr.row.size() != (n + 1))
//     {
//         printf("Error: sizes of row and c are incompatible\n");
//         exit(1);
//     }
//     else if (csr.col.size() != nz)
//     {
//         printf("Error: sizes of col and val are incompatible\n");
//         exit(1);
//     }
//     else if (((n * n) - 2 * n) < (2 * nz)) // n^2 > n + 2*nz,  resulting in CSR taking more space than dense matrix representation
//     {
//         printf("Error: CSR requires more space than dense matrix representation \n Use dense matrix implementation instead...\n");
//         exit(1);
//     }
//     else if ((csrM.row.size() != (n + 1)) || (csrM.col.size() != nz) || (csrM.val.size() != nz))
//     // in case the input graph is its minor and no further compression is possible, the dimensions of the compressed vectors should be equal to the dimensions of the input vectors
//     {
//         printf("Error: at least one of the compressed vectors doesn't have enough allocated space \n");
//         exit(1);
//     }

//     size_t nclus = 0;

//     numClusters(nclus, c); // find the number of distinct clusters

//     size_t end, localCount;
//     std::atomic<uint32_t> allCount(0);                        // common atomic variable to keep track of the index of the col and val CSR output
//     bool clusterHasElements = 0;                              // flag that denotes if this cluster has any elements
//     std::vector<std::atomic<uint32_t>> auxValueVector(nclus); // auxiliary vector that will contain all the non-zero values of each cluster (element of rowM)

//     csrM.row.resize(nclus + 1); // resize vector to the number of clusters

//     for (size_t id = 1; id < (nclus + 1); id++) // cluster ids start from 1
//     {
//         csrM.row[id - 1] = allCount.load();

//         cilk_for(size_t i = 0; i < nclus; i++)
//         {
//             auxValueVector[i] = 0; // reset auxiliary vector
//         }

//         clusterHasElements = 0;

//         cilk_for(size_t i = 0; i < n; i++)
//         {
//             if (id != c[i]) // c[i]: cluster of row i of colCompressed/row
//                 continue;

//             end = csr.row[i + 1];

//             if (csr.row[i] == end) // this row has no non-zero elements
//                 continue;          // OpenCilk allows work stealing, so continue is encouraged to make the thread finish earlier and steal work from other threads

//             if (!clusterHasElements) // declare that this row has non-zero elements only once
//                 clusterHasElements = 1;

//             for (size_t j = csr.row[i]; j < end; j++)
//             {
//                 auxValueVector[c[csr.col[j]] - 1] += (csr.val[j]); // compress cols by summing the values of each cluster to the column the cluster id points to
//             }
//         }

//         if (!clusterHasElements)
//             continue;

//         cilk_for(size_t i = 0; i < nclus; i++)
//         {
//             if (auxValueVector[i] == 0)
//                 continue;

//             localCount = allCount.fetch_add(1); // returns the value of allCount and then increments it by 1 atomically
//             csrM.val[localCount] = auxValueVector[i];
//             csrM.col[localCount] = i;
//         }
//     }

//     csrM.row[nclus] = allCount;
//     csrM.col.resize(allCount);
//     csrM.val.resize(allCount);
// }

void GMopenCilk(CSR &csrM, const CSR &csr, const std::vector<size_t> &c)
{
    size_t n = c.size();
    size_t nz = csr.val.size();

    // the row vector should be of size equal to the rows of the matrix, plus one (the latter is an implementation specific)
    // we assume that the c vector's size is set correctly equal to the number of nodes of the input grapg / rows of the input matrix
    if (csr.row.size() != (n + 1))
    {
        printf("Error: sizes of row and c are incompatible\n");
        exit(1);
    }
    else if (csr.col.size() != nz)
    {
        printf("Error: sizes of col and val are incompatible\n");
        exit(1);
    }
    else if (((n * n) - 2 * n) < (2 * nz)) // n^2 > n + 2*nz,  resulting in CSR taking more space than dense matrix representation
    {
        printf("Error: CSR requires more space than dense matrix representation \n Use dense matrix implementation instead...\n");
        exit(1);
    }
    else if ((csrM.row.size() != (n + 1)) || (csrM.col.size() != nz) || (csrM.val.size() != nz))
    // in case the input graph is its minor and no further compression is possible, the dimensions of the compressed vectors should be equal to the dimensions of the input vectors
    {
        printf("Error: at least one of the compressed vectors doesn't have enough allocated space \n");
        exit(1);
    }

    size_t nclus = 0;

    numClusters(nclus, c); // find the number of distinct clusters

    size_t end, localCount;
    std::atomic<uint32_t> allCount(0);                   // common atomic variable to keep track of the index of the col and val CSR output
    bool clusterHasElements = 0;                         // flag that denotes if this cluster has any elements
    std::vector<std::atomic<uint32_t>> commonAux(nclus); // common auxiliary vector to save compressed row and col of a single cluster

    csrM.row.resize(nclus + 1); // resize vector to the number of clusters

    for (size_t id = 1; id < (nclus + 1); id++) // cluster ids start from 1
    {
        csrM.row[id - 1] = allCount.load();

        cilk_for(size_t i = 0; i < nclus; i++)
        {
            commonAux[i] = 0; // reset auxiliary vector
        }

        clusterHasElements = 0;
        cilk_scope
        {
            std::vector<uint32_t> auxValueVector(nclus); // local auxiliary vector to save compressed row and col of a single cluster
            cilk_for(size_t i = 0; i < n; i++)
            {
                if (id != c[i]) // c[i]: cluster of row i of colCompressed/row
                    continue;

                end = csr.row[i + 1];

                if (csr.row[i] == end) // this row has no non-zero elements
                    continue;          // OpenCilk allows work stealing, so continue is encouraged to make the thread finish earlier and steal work from other threads

                if (!clusterHasElements) // declare that this row has non-zero elements only once
                    clusterHasElements = 1;

                for (size_t j = csr.row[i]; j < end; j++)
                {
                    auxValueVector[c[csr.col[j]] - 1] += (csr.val[j]); // compress cols by summing the values of each cluster to the column the cluster id points to
                }
            }

            if (clusterHasElements)
            {
                cilk_for(size_t i = 0; i < nclus; i++)
                {
                    if (auxValueVector[i] == 0)
                        continue;
                    commonAux[i].fetch_add(auxValueVector[i]); // combine the local auxiliary vectors in the common auxiliary vector
                                                               // produces less traffic than having the auxValueVector as shared and atomic
                }
            }
        }

        if (!clusterHasElements) // executed by a single thread (the main thread)
            continue;

        cilk_for(size_t i = 0; i < nclus; i++)
        {
            if (commonAux[i] == 0)
                continue;

            localCount = allCount.fetch_add(1); // returns the value of allCount and then increments it by 1 atomically
            csrM.val[localCount] = commonAux[i];
            csrM.col[localCount] = i;
        }
    }

    csrM.row[nclus] = allCount;
    csrM.col.resize(allCount);
    csrM.val.resize(allCount);
}
