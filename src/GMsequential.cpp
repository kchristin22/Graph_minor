#include "GMsequential.hpp"

inline void numClusters(size_t &nclus, const std::vector<size_t> &c)
{
    size_t n = c.size();

    size_t max = 0; // ids start from 1

    for (size_t i = 0; i < n; i++)
    {
        max = (c[i] > max) ? c[i] : max;
    }

    nclus = max;
}

inline void numClustersGen(size_t &nclus, const std::vector<size_t> &c)
{
    size_t n = c.size();
    std::vector<size_t> discreteClus(n, 0); // vector where the ith element is 1 if cluster i has any nodes

    for (size_t i = 0; i < n; i++)
    {
        discreteClus[c[i] - 1] = 1; // in dense representation, we assume that there is no ith row and column that are both zero so we know that all ids included in c exist in A
                                    // we can atomically add 1, instead, to the cluster of the ith row to know how many nodes are in each cluster
    }

    nclus = 0;
    for (size_t i = 0; i < n; i++)
    {
        if (discreteClus[i] == 0)
            continue;
        nclus++;
    }
}

void seq(std::vector<int> &M, const std::vector<int> &A, const std::vector<size_t> &c)
{
    if (A.size() != (c.size() * c.size())) // the input should be squared and the c vector should be equal to one of its dimensions
    {
        printf("Error: sizes of A and c are incompatible\n");
        exit(1);
    }
    size_t n = c.size();
    size_t nclus;

    numClusters(nclus, c); // get the number of clusters

    int x;
    std::vector<int> colCompressed(n * nclus, 0); // intermediate vector to store the compressed column version of the input matrix

    for (size_t i = 0; i < n; i++) // access A and colCompressed by row to take advantage of cache locality
    {
        for (size_t j = 0; j < n; j++)
        {
            x = A[i * n + j];
            if (x == 0)
                continue;
            colCompressed[i * nclus + c[j] - 1] += x; // compress cols by summing the values of each cluster to the column the cluster id points to
        }
    }

    // Note: if ids are continuous they don't need sorting to remove rows corresponding to missing clusters,
    // otherwise they do to avoid adding false nodes to the graph minor

    M.resize(nclus * nclus, 0); // now that nclus is knwon, we can resize the output to the correct size

    for (size_t id = 1; id < (nclus + 1); id++) // cluster ids start from 1
    {                                           // access one row of M at a time to take advantage of cache locality
        for (size_t i = 0; i < n; i++)
        {
            if (id != c[i]) // c[i]: cluster of row i of colCompressed/row
                continue;

            for (size_t j = 0; j < nclus; j++)
            {
                M[(id - 1) * nclus + j] += colCompressed[i * nclus + j]; // compress rows of colCompressed by summing the rows of each cluster to the row of M the cluster id points to
            }
        }
    }
};

// this is the main function used for sparse matrix implementation
void seq(CSR &csrM, const CSR &csr, const std::vector<size_t> &c)
{
    // the row vector should be of size equal to the rows of the matrix, plus one (the latter is an implementation specific)
    // we assume that the c vector's size is set correctly equal to the number of nodes of the input grapg / rows of the input matrix
    if (csr.row.size() != (c.size() + 1))
    {
        printf("Error: sizes of row and c are incompatible\n");
        exit(1);
    }
    else if (csr.col.size() != csr.val.size())
    {
        printf("Error: sizes of col and val are incompatible\n");
        exit(1);
    }
    else if (csr.row.size() == (csr.col.size() + 1)) // the number of non-zero elements is equal to the number of rows,
                                                     // resulting in CSR taking more space than dense matrix representation
    {
        printf("Error: CSR requires more space than dense matrix representation \n Use dense matrix implementation instead...\n");
        exit(1);
    }
    else if ((csrM.row.size() != csr.row.size()) || (csrM.col.size() != csr.col.size()) || (csrM.val.size() != csr.val.size()))
    // in case the input graph is its minor and no further compression is possible, the dimensions of the compressed vectors should be equal to the dimensions of the input vectors
    {
        printf("Error: at least one of the compressed vectors doesn't have enough allocated space \n");
        exit(1);
    }

    size_t n = c.size();
    size_t nclus;

    numClusters(nclus, c); // get the number of clusters

    size_t end, allCount = 0;                    // loop variables
    bool clusterHasElements = 0;                 // flag that denotes if the cluster has any elements
    std::vector<uint32_t> auxValueVector(nclus); // auxiliary vector that contains all the values, zero and non-zero, of a cluster
    csrM.row.resize(nclus + 1);                  // resize vector to the number of clusters, plus one (implementation specific)

    for (size_t id = 1; id < (nclus + 1); id++) // cluster ids start from 1
    {
        csrM.row[id - 1] = allCount;     // store offset of this rowM element
        auxValueVector.assign(nclus, 0); // reset auxiliary vector
        clusterHasElements = 0;          // reset flag

        for (size_t i = 0; i < n; i++)
        {
            if (id != c[i]) // c[i]: cluster of row i of colCompressed/row
                continue;

            end = csr.row[i + 1]; // if row.size()=n was true, we should had set end=nz in the final iteration

            if (csr.row[i] == end) // this row has no non-zero elements
                continue;
            if (!clusterHasElements) // declare that this row has non-zero elements, only one set of the flag needed
                clusterHasElements = 1;

            for (size_t j = csr.row[i]; j < end; j++)
            {
                auxValueVector[c[csr.col[j]] - 1] += csr.val[j]; // compress cols by summing the values of each cluster to the column the cluster id points to (-1, as ids start from 1)
            }
        }
        if (!clusterHasElements) // the cluster has no non-zero elements
            continue;
        for (size_t i = 0; i < nclus; i++) // store the connections of this cluster to all the rest
        // access the vectors sequentially to take advantage of cache locality
        {
            if (auxValueVector[i] == 0) // CSR stores only the non-zero elements
                continue;

            csrM.val[allCount] = auxValueVector[i]; // store value of the weighted edge
            csrM.col[allCount] = i;                 // store column/cluster id of the weighted edge
            allCount++;                             // update index
        }
    }
    csrM.row[nclus] = allCount; // last element of rowM is the number of non-zero elements of valM and colM
    csrM.col.resize(allCount);  // resize vectors to their actual size
    csrM.val.resize(allCount);
}

void seqDenseCSR(std::vector<size_t> &rowM, std::vector<size_t> &colM, std::vector<uint32_t> &valM,
                 const std::vector<size_t> &row, const std::vector<size_t> &col, const std::vector<uint32_t> &val, const std::vector<size_t> &c)
{
    // the row vector should be of size equal to the rows of the matrix, plus one (the latter is an implementation specific)
    // we assume that the c vector's size is set correctly equal to the number of nodes of the input grapg / rows of the input matrix
    if (row.size() != (c.size() + 1))
    {
        printf("Error: sizes of row and c are incompatible\n");
        exit(1);
    }
    else if (col.size() != val.size())
    {
        printf("Error: sizes of col and val are incompatible\n");
        exit(1);
    }
    else if (row.size() == col.size()) // the number of non-zero elements is equal to the number of rows,
                                       // resulting in CSR taking more space than dense matrix representation
    {
        printf("Error: CSR requires more space than dense matrix representation \n Use dense matrix implementation instead...\n");
        exit(1);
    }
    else if ((rowM.size() != row.size()) || (colM.size() != col.size()) || (valM.size() != val.size()))
    // in case the input graph is its minor and no further compression is possible, the dimensions of the compressed vectors should be equal to the dimensions of the input vectors
    {
        printf("Error: at least one of the compressed vectors doesn't have enough allocated space \n");
        exit(1);
    }

    size_t n = c.size();
    size_t nclus;

    numClusters(nclus, c); // get the number of clusters

    std::vector<int> colCompressed(n * nclus, 0); // intermediate vector to store the compressed column version of the input matrix

    for (size_t i = 0; i < n; i++) // access row, col, val and colCompressed sequentially to take advantage of cache locality
    {
        for (size_t j = row[i]; j < row[i + 1]; j++)
        {
            colCompressed[i * nclus + c[col[j]] - 1] += val[j]; // compress cols by summing the values of each cluster to the column the cluster id points to
        }
    }

    // Note: if ids are continuous they don't need sorting to remove rows corresponding to missing clusters,
    // otherwise they do to avoid adding false nodes to the graph minor

    std::vector<int> M(nclus * nclus, 0);

    for (size_t id = 1; id < (nclus + 1); id++) // cluster ids start from 1
    {                                           // access one row of M at a time to take advantage of cache locality
        for (size_t i = 0; i < n; i++)
        {
            if (id != c[i]) // c[i]: cluster of row i of colCompressed/row
                continue;

            for (size_t j = 0; j < nclus; j++)
            {
                M[(id - 1) * nclus + j] += colCompressed[i * nclus + j]; // compress rows of colCompressed by summing the rows of each cluster to the row of M the cluster id points to
            }
        }
    }

    size_t count = 0;
    int x;
    rowM.resize(nclus + 1);
    for (size_t i = 0; i < nclus; i++) // access rowM, colM, valM and M sequentially to take advantage of cache locality
    {
        rowM[i] = count;
        for (size_t j = 0; j < nclus; j++)
        {
            x = M[i * nclus + j];
            if (x == 0)
                continue;
            valM[count] = x;
            colM[count] = j;
            count++;
        }
    }
    rowM[nclus] = count;
    colM.resize(count);
    valM.resize(count);
}