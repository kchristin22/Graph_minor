# Graph_minor

## (!!!) There's an error in the benchmarks in the report: instead of pthreads, OpenMP was run. The fixed tests are below.

### **About:**

The purpose of this project is to calculate the graph minor of a given graph, based on an input vector denoting the clusters each original node belongs to. This vector contains the ids of the
clusters, ranging from 1 to the amount of discrete clusters, in the equivalent positions of the corresponding nodes. No assumptions on the nature of the input graph are made. However, since most of the big networks are sparse, the primary emphasis lies on sparse matrix representations, and partic-ularly in the Compressed Sparse Row (CSR) format. It is important to note that this report does not address any preprocessing steps required to convert the input adjacency matrix into CSR format.

### **Benchmarks:**

Comparison of implementations using the median time in us:

| Data set | Matlab |  Sequential | OpenMP (4 threads) | Pthreads  (4 threads) | OpenCilk  (4 threads)|
|---------|-----------|-----------|----------|---------| -----------|
| ri2010  | 4500 | 34289 | 9767 | 33618 |  19133 |
| nj2010  | 56600 |1584908 | 339045 | 545741 | 544979 |
| pa2010  | 110000 |6037295 | 2060738 | 2739690 | 2103721 |
| ca2010  | 187400 |15898973  | 5423717 | 6440966 | 4838647 | 
| tx2010 | 231300 |25460247 | 8036885 | 9773966 | 7118827 |

OpenCilk manages to outperform OpenMP in larger data sets, despite the use of a common atomic vector. This is assumed to stem from its work-stealing properties.

![image](https://github.com/kchristin22/Graph_minor/assets/74819775/58646b43-75ec-4073-94b8-110b17959da5)

 The `us2010` dataset was never managed to be parsed correctly or quicky enough.
       

