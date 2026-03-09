#ifndef MM_HPP
#define MM_HPP

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "readgri.hpp"
#include "quad.hpp"

/* 2x2 Jacobian matrix for the linear mapping 
 * from the reference triangle to physical element k.
 */
Eigen::Matrix2d Jacobian(const GriMesh& mesh, int elem);

/* Computes the Mass Matrix on the reference triangle:
 * M_ref_ij = integral( phi_i * phi_j ) over reference element.
 */
Eigen::MatrixXd computeRefMassMatrix(int order);

/* Assembles the global Block-Diagonal Mass Matrix.
 * For DG, this is a sparse matrix where each block is (Np x Np).
 * Layout corresponds to [element * Np + node].
 */
Eigen::SparseMatrix<double> computeGlobalMassMatrix(const GriMesh& mesh, int order);

#endif
