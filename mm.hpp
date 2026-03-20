#ifndef MM_HPP
#define MM_HPP

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "readgri.hpp"
#include "quad.hpp"

/* 2x2 Jacobian matrix for the affine mapping built from element corner nodes.
 * Kept for legacy/debug usage; curved geometry uses geometry.hpp utilities.
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

/**
 * Apply the block-diagonal mass matrix inverse in-place to a DG residual vector R.
 *
 * The DG time equation is  M * dU/dt = -R  (R = raw residual from calcRes).
 * The correct time update is therefore  U_new = U - dt * M^{-1} * R.
 * This function overwrites R with  M^{-1} * R  so that the caller can write
 *   U -= dt * R  directly.
 *
 * For curved elements, M_k is assembled with quadrature using geometric Jacobian
 * at each quadrature point, and M_k^{-1} is cached per element.
 *
 * R layout: [var * Ne * Np + elem * Np + basis]
 */
void applyInverseMassMatrix(const GriMesh& mesh, double* R, int order);

#endif
