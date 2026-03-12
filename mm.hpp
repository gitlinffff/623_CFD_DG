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

/**
 * Apply the block-diagonal mass matrix inverse in-place to a DG residual vector R.
 *
 * The DG time equation is  M * dU/dt = -R  (R = raw residual from calcRes).
 * The correct time update is therefore  U_new = U - dt * M^{-1} * R.
 * This function overwrites R with  M^{-1} * R  so that the caller can write
 *   U -= dt * R  directly.
 *
 * For a linear (affine) element:  M_k = det(J_k) * M_ref
 * so  M_k^{-1} = (1/det(J_k)) * M_ref^{-1}.
 * M_ref^{-1} is the same for every element and is cached after the first call.
 *
 * R layout: [var * Ne * Np + elem * Np + basis]
 */
void applyInverseMassMatrix(const GriMesh& mesh, double* R, int order);

#endif
