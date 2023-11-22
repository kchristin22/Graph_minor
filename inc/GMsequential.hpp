#pragma once

#include "coo_to_csr.hpp"

/* Calculates the number of discrete clusters, in the case of having contiguous ids
 * @params:
 *   nclus (output): number of discrete clusters
 *   c (input): vector containing the cluster of each node of the input graph
 */
inline void numClusters(size_t &nclus, const std::vector<size_t> &c);

/* Calculates the number of discrete clusters. This function is an idea of what counting the ids would look like if they weren't contiguous.
 * @params:
 *   nclus (output): number of discrete clusters
 *   c (input): vector containing the cluster of each node of the input graph
 */
inline void numClustersGen(size_t &nclus, const std::vector<size_t> &c);

/* Dense matrix implementation
 * @params:
 *   M (output): dense matrix representation of the adjacency matrix corresponding to the graph minor of the input graph
 *   A (input): dense matrix representation of the adjacency matrix of the input graph
 *   c (input): vector containing the cluster of each node of the input graph
 */
void seq(std::vector<int> &M, const std::vector<int> &A, const std::vector<size_t> &c);

/* CSR matrix implementation, this is the function used
 * @params:
 *   csrM (output): csr representation of the adjacency matrix corresponding to the graph minor of the input graph
 *   csr (input): csr representation of the adjacency matrix of the input graph
 *   c (input): vector containing the cluster of each node of the input graph
 */
void seq(CSR &csrM, const CSR &csr, const std::vector<size_t> &c);

/* Calculates CSR result using dense matrix implementation in intermediate steps
 * @params:
 *   rowM (output): CSR row vector of the adjacency matrix corresponding to the graph minor of the input graph
 *   colM (output): CSR col vector of the adjacency matrix corresponding to the graph minor of the input graph
 *   valM (output): CSR val vector of the adjacency matrix corresponding to the graph minor of the input graph
 *   row (input): CSR row vector of the adjacency matrix of the input graph
 *   col (input): CSR col vector of the adjacency matrix of the input graph
 *   val (input): CSR val vector of the adjacency matrix of the input graph
 *   c (input): vector containing the cluster of each node of the input graph
 */
void seqDenseCSR(std::vector<size_t> &rowM, std::vector<size_t> &colM, std::vector<uint32_t> &valM,
                 const std::vector<size_t> &row, const std::vector<size_t> &col, const std::vector<uint32_t> &val, const std::vector<size_t> &c);
