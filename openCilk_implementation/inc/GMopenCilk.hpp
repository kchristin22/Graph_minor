#pragma once

#include <cilk/cilk.h>
#include "coo_to_csr.hpp"

/* Calculates the number of discrete clusters, in the case of having contiguous ids
 * @params:
 *   nclus (output): number of discrete clusters
 *   c (input): vector containing the cluster of each node of the input graph
 */
inline void numClusters(size_t &nclus, const std::vector<size_t> &c);

/* The sequential algorithm implemented using OpenCilk
 * @params:
 *   csrM (output): CSR representation of the adjacency matrix corresponding to the graph minor of the input graph
 *   csr (input): CSR representation of the adjacency matrix of the input graph
 *   c (input): vector containing the cluster of each node of the input graph
 */
void GMopenCilk(CSR &csrM, const CSR &csr, const std::vector<size_t> &c);
