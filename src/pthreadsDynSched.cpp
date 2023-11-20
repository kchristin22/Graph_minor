#include "GMpthreads.hpp"

bool finishedTasks = 0;

void *producer(void *args)
{
    Task task = *(Task *)args;

    // add check of max queue size

    for (size_t i = task.chunkStart; i < task.chunkEnd; i += ELEMENTS_PER_CACHE_LINE) // each producer thread takes a chunk start and end and puts a cacheline in
    {
        task.nclusArgs.start = i;
        task.nclusArgs.end = i + ELEMENTS_PER_CACHE_LINE;
        pthread_mutex_lock(task.queueMutex);
        task.queue->push(task);
        pthread_mutex_unlock(task.queueMutex);
        pthread_cond_signal(task.queueNotEmpty);
    }

    return nullptr;
}

void *consumer(void *args)
{
    Task *task = (Task *)args;

    while (1)
    {

        pthread_mutex_lock(task->queueMutex);
        if (task->queue->empty())
        {
            if (finishedTasks)
            {
                pthread_mutex_unlock(task->queueMutex);
                return nullptr;
            }
            pthread_cond_wait(task->queueNotEmpty, task->queueMutex);
        }

        Task tmpTask = task->queue->front();
        task->queue->pop();
        pthread_mutex_unlock(task->queueMutex);

        tmpTask.fn(&tmpTask.nclusArgs);
    }

    return nullptr;
}

inline void *fnNumClusters(void *args)
{
    nclusThread nclusArgs = *(nclusThread *)args;
    size_t n = nclusArgs.c.size();
    std::vector<size_t> discreetClus(n, 0); // vector where the ith element is a if cluster i has a nodes

    for (size_t i = nclusArgs.start; i < nclusArgs.end; i++)
    {
        discreetClus[(nclusArgs.c[i] - 1)] = 1; // we assume that there is no ith row and column that are both zero so we know that all ids included in c exist in A
                                                // we can atomically add 1, instead, to the cluster of the ith row to know how many nodes are in each cluster
    }

    for (size_t i = nclusArgs.start; i < nclusArgs.end; i++)
    {
        if (discreetClus[i] == 0)
            continue;
        nclusArgs.nclus.fetch_add(1);
    }

    return nullptr;
}

void GMpthreads(CSR &csrM, const CSR &csr, const std::vector<size_t> &c)
{
    if (csr.row.size() != c.size())
    {
        printf("Error: sizes of row and c are incompatible\n");
        exit(1);
    }
    else if (csr.col.size() != csr.val.size())
    {
        printf("Error: sizes of col and val are incompatible\n");
        exit(1);
    }
    else if (csr.row.size() == csr.col.size())
    {
        printf("Error: CSR requires more space than dense matrix representation \n Use dense matrix implementation instead...\n");
        exit(1);
    }
    else if ((csrM.row.size() != csr.row.size()) || (csrM.col.size() != csr.col.size()) || (csrM.val.size() != csr.val.size()))
    // if the input graph is its the minor, the dimensions of the compressed vectors should be equal to the dimensions of the input vectors
    {
        printf("Error: at least one of the compressed vectors doesn't have enough allocated space \n");
        exit(1);
    }

    size_t n = c.size();
    size_t nz = csr.val.size();

    size_t cacheLines = n / 16; // how many cache lines the array fills:
                                // a chunk of x elements of an int array (4*x bytes) will be equal to a cache line (64 bytes) to avoid false sharing
    size_t numThreads = 4;      // consider edge case where we need + 1 thread

    if (!cacheLines)
    {
        numThreads = 1;
        cacheLines = 1;
    }
    else if (cacheLines <= numThreads) // the array of the smallest type (int < size_t) fits in <=4 cache lines
    {
        numThreads = cacheLines;
    }

    size_t chunk = cacheLines * ELEMENTS_PER_CACHE_LINE / numThreads;
    size_t lastThreadChunk = chunk;

    if (n > (cacheLines * ELEMENTS_PER_CACHE_LINE)) // n does not fill in an integer number of cachelines
    {
        lastThreadChunk = chunk + (n - cacheLines * ELEMENTS_PER_CACHE_LINE); // the last thread gets its chunk + the remaining items
    }
    else // the array fits in a single cache line
        lastThreadChunk = n;

    std::atomic<size_t> nclus{0}; // if we used different variables for each thread and then calculate their sum,
                                  // since this vector of variables would fit in a single cache line,
                                  // the whole vector would be updated each time an element is changed, so it would be as or more costly

    std::vector<nclusThread> nclusArgs;
    nclusArgs.reserve(numThreads);

    // task queue

    printf("Num threads: %ld\n", numThreads);

    std::queue<Task> queue;
    pthread_mutex_t queueMutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_cond_t queueNotEmpty = PTHREAD_COND_INITIALIZER;

    std::vector<Task> tasks;
    tasks.reserve(numThreads);

    pthread_t nclusThreads[numThreads], cons[numThreads];

    for (size_t i = 0; i < numThreads; i++)
    {
        size_t start = i * chunk;
        size_t end = (i == (numThreads - 1)) ? (start + lastThreadChunk) : (start + chunk);
        nclusArgs.push_back({0, 0, c, nclus}); // nclusArgs used to include an atomic size_t for the number of clusters
        tasks.push_back({&queue, &queueMutex, &queueNotEmpty, start, end, nclusArgs[i], fnNumClusters});
        pthread_create(&nclusThreads[i], NULL, producer, (void *)&tasks[i]);
    }

    consThread consTasks[numThreads] = {{&queue,
                                         &queueMutex,
                                         &queueNotEmpty},
                                        {&queue,
                                         &queueMutex,
                                         &queueNotEmpty},
                                        {&queue,
                                         &queueMutex,
                                         &queueNotEmpty},
                                        {&queue,
                                         &queueMutex,
                                         &queueNotEmpty}};

    for (size_t i = 0; i < numThreads; i++)
    {
        pthread_create(&cons[i], NULL, consumer, (void *)&consTasks[i]);
    }

    for (pthread_t i : nclusThreads)
    {
        pthread_join(i, NULL);
    }

    finishedTasks = 1;

    for (pthread_t i : cons)
    {
        pthread_cond_broadcast(&queueNotEmpty); // broadcasting signals is used to wake any consumer threads that are waiting
        pthread_join(i, NULL);
    }

    printf("Number of clusters: %ld\n", nclus.load());
    printf("I'm in the dynamic\n");

    return;
}