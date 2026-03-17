#include "mm.hpp"
#include <cmath>
#include <stdexcept>

extern "C" {
	void shapeL(double *xref, int p, double **pphi);
}

Eigen::Matrix2d Jacobian(const GriMesh& mesh, int elem) {

	Eigen::Matrix2d J;

	int v1_idx = mesh.E[elem * 3 + 0];
	int v2_idx = mesh.E[elem * 3 + 1];
	int v3_idx = mesh.E[elem * 3 + 2];

	double x1 = mesh.V[v1_idx * 2 + 0], y1 = mesh.V[v1_idx * 2 + 1];
	double x2 = mesh.V[v2_idx * 2 + 0], y2 = mesh.V[v2_idx * 2 + 1];
	double x3 = mesh.V[v3_idx * 2 + 0], y3 = mesh.V[v3_idx * 2 + 1];

	J << x2-x1, x3-x1,
	     y2-y1, y3-y1;
	return J;
}

Eigen::MatrixXd computeRefMassMatrix(int order) {
	int Np = (order + 1) * (order + 2) / 2;
	Eigen::MatrixXd M_ref = Eigen::MatrixXd::Zero(Np, Np);

	QuadratureRule quad = getQuadratureRule(order);

	double* phi = nullptr;
	double xref[2];

	for (int q = 0; q < quad.nq; ++q) {
		xref[0] = quad.xq[2 * q];
		xref[1] = quad.xq[2 * q + 1];

		shapeL(xref, order, &phi);

		for (int i = 0; i < Np; ++i)
			for (int j = 0; j < Np; ++j) {
				M_ref(i, j) += quad.wq[q] * phi[i] * phi[j];
			}
	}

	free(phi);
	return M_ref;
}

Eigen::SparseMatrix<double> computeGlobalMassMatrix(const GriMesh& mesh, int order) {
	int Np = (order + 1) * (order + 2) / 2;
	int totalSize = mesh.Ne * Np;

	Eigen::SparseMatrix<double> M(totalSize, totalSize);
	Eigen::MatrixXd M_ref = computeRefMassMatrix(order);

	std::vector<Eigen::Triplet<double>> tripletList;
	tripletList.reserve(mesh.Ne * Np * Np);

	for (int k = 0; k < mesh.Ne; ++k) {
		Eigen::Matrix2d J = Jacobian(mesh, k);
		double detJ = J.determinant();
		if (detJ <= 1e-15) {
			throw std::runtime_error("Negative or zero detJ at element " + std::to_string(k));
		}

		int blockStart = k * Np;
		for (int i = 0; i < Np; ++i) {
			for (int j = 0; j < Np; ++j) {
				double val = detJ * M_ref(i, j);
				tripletList.push_back(Eigen::Triplet<double>(blockStart + i, blockStart + j, val));
			}
		}
	}

	M.setFromTriplets(tripletList.begin(), tripletList.end());
	M.makeCompressed();
	return M;
}

void applyInverseMassMatrix(const GriMesh& mesh, double* R, int order) {
	int Np = (order + 1) * (order + 2) / 2;

	static Eigen::MatrixXd M_ref_inv_cache;
	static int cached_order = -1;
	if (cached_order != order) {
		M_ref_inv_cache = computeRefMassMatrix(order).inverse();
		cached_order = order;
	}
	const Eigen::MatrixXd& Minv = M_ref_inv_cache;

	#pragma omp parallel for schedule(static)
	for (int k = 0; k < mesh.Ne; ++k) {
		Eigen::Matrix2d J = Jacobian(mesh, k);
		double detJ = std::abs(J.determinant());
		if (detJ < 1e-30) detJ = 1e-30;
		double inv_detJ = 1.0 / detJ;

		Eigen::VectorXd R_elem(Np);
		for (int var = 0; var < 4; ++var) {
			int base = (var * mesh.Ne + k) * Np;
			for (int i = 0; i < Np; ++i)
				R_elem(i) = R[base + i];
			Eigen::VectorXd R_new = inv_detJ * Minv * R_elem;
			for (int i = 0; i < Np; ++i)
				R[base + i] = R_new(i);
		}
	}
}
