#ifndef MM_HPP
#define MM_HPP

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "readgri.hpp"
#include "quad.hpp"

Eigen::Matrix2d Jacobian(const GriMesh& mesh, int elem);

Eigen::MatrixXd computeRefMassMatrix(int order);

Eigen::SparseMatrix<double> computeGlobalMassMatrix(const GriMesh& mesh, int order);

void applyInverseMassMatrix(const GriMesh& mesh, double* R, int order);

#endif
