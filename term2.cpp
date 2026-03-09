#include <vector>
#include <Eigen/Dense>
#include "physics.hpp"
#include "quad.hpp"
#include "readgri.hpp"
#include "problem.hpp"

//  check and test this function
void addTerm2(const GriMesh& mesh, std::vector<double>& R, int order,
	            const std::vector<BasisEval>& phiq, const double* U,
							const ProblemParams& params) {
	int Np = (order + 1) * (order + 2) / 2;
	QuadratureRule quad = getQuadratureRule(order);

	for (int k = 0; k < mesh.Ne; ++k) {
		Eigen::Matrix2d J = Jacobian(mesh, k);
		double detJ = J.determinant();
		Eigen::Matrix2d J_inv = J.inverse();

		for (int q = 0; q < quad.nq; ++q) {
			// Interpolate the 4 variables at this quad point
			Eigen::Vector4d uq = Eigen::Vector4d::Zero();
			for (int i = 0; i < Np; ++i) {
				double phi_val = phiq[i * quad.nq + q].phi;
				for (int var = 0; var < 4; ++var) {
					// U  [variable][element][node]
					uq(var) += U(var*mesh.Ne*Np+k*Np+i) * phi_val;
				}
			}

			// Evaluate Physical Flux at this quad point: [Fx, Fy] each is a Vector4d
			Eigen::Matrix<double, 4, 2> phyF = vec_phyFlux(uq.data(), params.gammad)

			// Loop over test functions to update the residual
			for (int i = 0; i < Np; ++i) {
				// Get gradient of phi
				Eigen::RowVector2d gradPhi_ref;
				gradPhi_ref << phiq[i * quad.nq + q].dphi_dxi,
											 phiq[i * quad.nq + q].dphi_deta;
				Eigen::RowVector2d gradPhi_phys = gradPhi_ref * J_inv;
				// loop over 4 variables
				for (int var = 0; var < 4; ++var) {
					double integral = (gradPhi_phys(0) * phyF(var, 0) + gradPhi_phys(1) * phyF(var, 1));
					R[(var*mesh.Ne+k)*Np+i] -= integral * detJ * quad.wq[q];
				}
			}
		}
	}
}
