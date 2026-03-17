#include "term2.hpp"
#include "physics.hpp"
#include "mm.hpp"

void addTerm2(const GriMesh& mesh, double* R, int order,
	            const std::vector<BasisEval>& phiq, const double* U,
							const ProblemParams& params) {
	int Np = (order + 1) * (order + 2) / 2;
	QuadratureRule quad = getQuadratureRule(order);

	#pragma omp parallel for schedule(static)
	for (int k = 0; k < mesh.Ne; ++k) {
		Eigen::Matrix2d J = Jacobian(mesh, k);
		double detJ = abs(J.determinant());
		Eigen::Matrix2d J_inv = J.inverse();

		for (int q = 0; q < quad.nq; ++q) {
			Eigen::Vector4d uq = Eigen::Vector4d::Zero();
			for (int i = 0; i < Np; ++i) {
				double phi_val = phiq[i * quad.nq + q].phi;
				for (int var = 0; var < 4; ++var) {
					uq(var) += U[var*mesh.Ne*Np+k*Np+i] * phi_val;
				}
			}

			Eigen::Matrix<double, 4, 2> phyF = vec_phyFlux(uq.data(), params.gammad);

			for (int i = 0; i < Np; ++i) {
				Eigen::RowVector2d gradPhi_ref;
				gradPhi_ref << phiq[i * quad.nq + q].dphi_dxi,
											 phiq[i * quad.nq + q].dphi_deta;
				Eigen::RowVector2d gradPhi_phys = gradPhi_ref * J_inv;
				for (int var = 0; var < 4; ++var) {
					double integral = (gradPhi_phys(0) * phyF(var, 0) + gradPhi_phys(1) * phyF(var, 1));
					R[(var*mesh.Ne+k)*Np+i] -= integral * detJ * quad.wq[q];
				}
			}
		}
	}
}
