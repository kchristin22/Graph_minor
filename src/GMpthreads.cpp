#include "GMpthreads.hpp"

pthread_barrier_t barrier; // barrier to synchronize threads, can also be passed as an argument to the threads

inline void calChunk(size_t &chunk, size_t &lastThreadChunk, const size_t n, const size_t min_chunk, const uint32_t numThreads)
{
    size_t cacheLines = n / min_chunk; // we assume that the minimum chunk given fills at least a whole cache line
                                       // f.i. a min chunk equal to 8 elements of a size_t array (8*8 bytes) will be equal to a cache line (64 bytes)

    if (!cacheLines) // the array fits in a single cache line
    {
        cacheLines = 1;
    }

    chunk = cacheLines * min_chunk / numThreads; // divide cache lines into chunks of equal size
    if (!chunk)                                  // if we have too many threads, then threads must share a cache line
    {
        chunk = n / numThreads; // we assign equal number of elements to each thread
    }

    lastThreadChunk = (n > (chunk * numThreads)) ? (chunk + (n - (chunk * numThreads))) : chunk; // if n does not fill in an integer number of cachelines,
                                                                                                 // the last thread gets its chunk + the remaining items
}

void *fnNumClusters(void *args)
{
    nclusThread *nclusArgs = (nclusThread *)args;

    size_t max = 0;

    for (size_t i = nclusArgs->start; i < nclusArgs->end; i++)
    {
        max = (nclusArgs->c[i] > max) ? nclusArgs->c[i] : max; // if this cluster id is greater than the current max, update the max
    }
    nclusArgs->nclus = max; // maximum id measured by the thread in that range of the vector c

    return nullptr;
}

void *fnThread(void *args)
{

    threadArgs *Args = (threadArgs *)args;
    size_t start = Args->start;
    size_t end = Args->end;

    if (end > (Args->c.size())) // check for an error in assigning the end variable (exceeding the size of the vector)
        end = Args->c.size();

    size_t startClus = Args->startClus; // range for vector of size nclus
    size_t endClus = Args->endClus;
    if (endClus > Args->nclus)
        endClus = Args->nclus;

    std::vector<uint32_t> auxValueVector(Args->nclus, 0); // local auxiliary vector to save compressed row and col of a single cluster

    size_t localCount = 0; // local storage of the index of the col and val CSR output

    for (size_t id = 1; id < (Args->nclus + 1); id++) // cluster ids start from 1
    {

        pthread_barrier_wait(&barrier); // synch threads before assigning the offset of the row vector

        Args->csrM.row[id - 1] = Args->allCount.load();

        auxValueVector.assign(Args->nclus, 0); // clear the vector for the next cluster

        for (size_t i = start; i < end; i++)
        {
            if (id != Args->c[i]) // c[i] : cluster of row i of colCompressed / row
                continue;

            for (size_t j = Args->csr.row[i]; j < Args->csr.row[i + 1]; j++)
            {
                auxValueVector[Args->c[Args->csr.col[j]] - 1] += (Args->csr.val[j]); // compress cols by summing the values of each cluster to the column the cluster id points to
            }
        }

        for (size_t i = 0; i < Args->nclus; i++)
        {
            if (auxValueVector[i] == 0)
                continue;
            Args->commonAux[i].fetch_add(auxValueVector[i]); // combine the local auxiliary vectors in the common auxiliary vector
                                                             // produces less traffic than having the auxValueVector as shared and atomic
        }

        pthread_barrier_wait(&barrier); // sync threads to have the commonAux vector filled correctly before using it for the col and val vectors
                                        // even if a thread didn't find any elements and its auxValueVector contains only zeros, it too must reach this point,
                                        // due to the barrier awaiting a specific number of threads to reach it
                                        // in addition, each thread has to execute its share of work(range) for the vectors of size nclus

        for (size_t i = startClus; i < endClus; i++)
        {
            if (Args->commonAux[i].load() == 0)
                continue;

            localCount = Args->allCount.fetch_add(1); // returns the value of allCount and then increments it by 1 atomically

            Args->csrM.val[localCount] = Args->commonAux[i];
            Args->csrM.col[localCount] = i;
            Args->commonAux[i] = 0; // clear the vector for the next cluster
        }
    }

    return nullptr;
}

void GMpthreads(CSR &csrM, const CSR &csr, const std::vector<size_t> &c, const uint32_t numThreads)
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

    size_t chunk, lastThreadChunk;
    calChunk(chunk, lastThreadChunk, n, ELEMENTS_PER_CACHE_LINE_SIZE_T, numThreads); // how many cache lines the array fills to avoid false sharing:
                                                                                     // c vector elements are size_t -> min chunk = cache line / sizeof(size_t)

    std::vector<size_t> nclusVector(numThreads, 0); // store the local max id measured by each thread
    std::vector<nclusThread> nclusArgs;             // thread arguments for calculating the number of clusters
    nclusArgs.reserve(numThreads);

    pthread_t nclusThreads[numThreads];

    for (size_t i = 0; i < numThreads; i++)
    {
        size_t start = i * chunk;
        size_t end = (i == (numThreads - 1)) ? (start + lastThreadChunk) : (start + chunk); // if it's the last thread, assign the last chunk
        nclusArgs.push_back({start, end, c, nclusVector[i]});
        pthread_create(&nclusThreads[i], NULL, fnNumClusters, (void *)&nclusArgs[i]);
    }

    for (pthread_t &i : nclusThreads)
    {
        pthread_join(i, NULL);
    }

    size_t nclus = 0; // initialize the number of clusters to 0 (the cluster id begin from 1)

    for (size_t i = 0; i < numThreads; i++)
    {
        nclus = (nclusVector[i] > nclus) ? nclusVector[i] : nclus; // if this cluster id is greater than the current max, update the max
    }

    // re-calculate cachelines and chunk sizes to cover for the INT vectors in that code block (uint32_t < size_t, so more of its elements fit in a cache line)
    calChunk(chunk, lastThreadChunk, n, ELEMENTS_PER_CACHE_LINE_INT, numThreads);

    size_t chunkClus, lastThreadChunkClus;
    calChunk(chunkClus, lastThreadChunkClus, nclus, ELEMENTS_PER_CACHE_LINE_INT, numThreads);

    std::atomic<size_t> allCount(0); // common atomic variable to keep track of the index of the col and val CSR output
    csrM.row.resize(nclus + 1);      // resize vector to the number of clusters

    pthread_t threads[numThreads];

    std::vector<threadArgs> argsVector;
    argsVector.reserve(numThreads);

    std::vector<std::atomic<uint32_t>> commonAux(nclus); // auxiliary vector that will contain all the non-zero values of each cluster (element of rowM)
    // the above part of calculating the nclus could not be moved inside the fnThread function,
    // as std::atomic objects are not copyable and thus a vector of them cannot be resized

    size_t start, startClus, end, endClus;

    pthread_barrier_init(&barrier, NULL, numThreads); // initialize barrier used to synchronize threads

    for (size_t i = 0; i < numThreads; i++)
    {
        start = i * chunk;
        end = (i == (numThreads - 1)) ? (start + lastThreadChunk) : (start + chunk);
        startClus = i * chunkClus;
        endClus = (i == (numThreads - 1)) ? (startClus + lastThreadChunkClus) : (startClus + chunkClus);
        argsVector.push_back({nclus, start, end, startClus, endClus, csr, c, commonAux, allCount, csrM});
        pthread_create(&threads[i], NULL, fnThread, (void *)&argsVector[i]);
    }

    for (pthread_t &i : threads)
    {
        pthread_join(i, NULL);
    }

    csrM.row[nclus] = allCount; // last element of rowM is the number of non-zero elements of valM and colM
    csrM.col.resize(allCount);
    csrM.val.resize(allCount);

    pthread_barrier_destroy(&barrier);
}