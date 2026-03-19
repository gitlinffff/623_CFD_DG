#include "term2.hpp"
#include "physics.hpp"
#include "mm.hpp"


struct ElementMetrics {
    std::vector<Eigen::Matrix2d> G; // store J_inv * detJ
    std::vector<double> detJ;
};


void precomputeCurvedMetrics(const GriMesh& mesh, int order,
	                          const std::vector<BasisEval>& phiq,
														std::map<int, ElementMetrics>& curved_metrics) {
	QuadratureRule quad = getQuadratureRule(order);

	for (int ib = 0; ib < mesh.num_boundary_faces; ++ib) {
		int bgroup = mesh.B2E[3 * ib + 2];
		if (classify_boundary_name(mesh.Bname[bgroup-1]) != DGBcType::WALL) continue;

		int elem = mesh.B2E[3 * ib + 0];
		int face = mesh.B2E[3 * ib + 1];

		// 1. Get Corner Indices from E
		int c0 = mesh.E[3 * elem + 0];
		int c1 = mesh.E[3 * elem + 1];
		int c2 = mesh.E[3 * elem + 2];

		double x_nodes[6], y_nodes[6];

		// 2. Map Corners (phi[0], phi[1], phi[2])
		// phi[0] @ (0,0), phi[1] @ (1,0), phi[2] @ (0,1)
		int corners[3] = {c0, c1, c2};
		for(int i = 0; i < 3; ++i) {
			x_nodes[i] = mesh.V[corners[i] * 2 + 0];
			y_nodes[i] = mesh.V[corners[i] * 2 + 1];
		}

		// 3. Map Mid-nodes (phi[3], phi[4], phi[5])
		int bstart = mesh.BedgeNodeOffset[ib];
		int mid_curved_idx = mesh.BedgeNodes[bstart + 2];
		double mx_c = mesh.V[mid_curved_idx * 2 + 0];
		double my_c = mesh.V[mid_curved_idx * 2 + 1];

		// Face 0: c0 -> c1 (Mid is phi[5])
		// Face 1: c1 -> c2 (Mid is phi[3])
		// Face 2: c2 -> c0 (Mid is phi[4])

		// Default all mid-nodes to linear averages first
		x_nodes[5] = 0.5 * (x_nodes[0] + x_nodes[1]); y_nodes[5] = 0.5 * (y_nodes[0] + y_nodes[1]);
		x_nodes[3] = 0.5 * (x_nodes[1] + x_nodes[2]); y_nodes[3] = 0.5 * (y_nodes[1] + y_nodes[2]);
		x_nodes[4] = 0.5 * (x_nodes[2] + x_nodes[0]); y_nodes[4] = 0.5 * (y_nodes[2] + y_nodes[0]);

		// Overwrite the ONE curved edge with the actual boundary node
		if (face == 0) { x_nodes[5] = mx_c; y_nodes[5] = my_c; }
		else if (face == 1) { x_nodes[3] = mx_c; y_nodes[3] = my_c; }
		else if (face == 2) { x_nodes[4] = mx_c; y_nodes[4] = my_c; }

		// 3. Compute Jacobian at each volume quadrature point
		ElementMetrics em;
		em.G.resize(quad.nq);
		em.detJ.resize(quad.nq);

		for (int q = 0; q < quad.nq; ++q) {
			double xi  = quad.xq[2*q];
			double eta = quad.xq[2*q+1];

			double dx_dxi = 0, dx_deta = 0, dy_dxi = 0, dy_deta = 0;
			for (int j = 0; j < 6; ++j) {
				dx_dxi  += x_nodes[j] * phiq[j * quad.nq + q].dphi_dxi;
				dx_deta += x_nodes[j] * phiq[j * quad.nq + q].dphi_deta;
				dy_dxi  += y_nodes[j] * phiq[j * quad.nq + q].dphi_dxi;
				dy_deta += y_nodes[j] * phiq[j * quad.nq + q].dphi_deta;
			}

			double detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;
			em.detJ[q] = std::abs(detJ);

			// Compute adjugate matrix (J_inv * detJ)
			Eigen::Matrix2d adj;
			adj <<  dy_deta, -dx_deta, 
						 -dy_dxi,  dx_dxi;
			em.G[q] = adj;
		}
		curved_metrics[elem] = em;
	}
}

void addTerm2(const GriMesh& mesh, double* R, int order,
	            const std::vector<BasisEval>& phiq, const double* U,
							const ProblemParams& params,
							const std::map<int, ElementMetrics>& curved_metrics) {
	int Np = (order + 1) * (order + 2) / 2;
	QuadratureRule quad = getQuadratureRule(order);

	#pragma omp parallel for schedule(static)
	for (int k = 0; k < mesh.Ne; ++k) {
		auto it = curved_metrics.find(k);
		bool is_curved = (it != curved_metrics.end());
		
		Eigen::Matrix2d J_const, G_const;
		if (!is_curved) {
			J_const = Jacobian(mesh, k);
			G_const = J.inverse() * abs(J.determinant());
		}

		for (int q = 0; q < quad.nq; ++q) {
			Eigen::Vector4d uq = Eigen::Vector4d::Zero();
			for (int i = 0; i < Np; ++i) {
				double phi_val = phiq[i * quad.nq + q].phi;
				for (int var = 0; var < 4; ++var) {
					uq(var) += U[var*mesh.Ne*Np+k*Np+i] * phi_val;
				}
			}

			Eigen::Matrix<double, 4, 2> phyF = vec_phyFlux(uq.data(), params.gammad);

			Eigen::Matrix2d current_G;
			if (is_curved) {
				current_G    = it->second.G[q]; 
			} else {
				current_G    = G_const;
			}

			for (int i = 0; i < Np; ++i) {
				Eigen::RowVector2d gradPhi_ref;
				gradPhi_ref << phiq[i * quad.nq + q].dphi_dxi,
											 phiq[i * quad.nq + q].dphi_deta;
				Eigen::RowVector2d gradPhi_phys = gradPhi_ref * current_G;
				for (int var = 0; var < 4; ++var) {
					double integral = (gradPhi_phys(0) * phyF(var, 0) + gradPhi_phys(1) * phyF(var, 1));
					R[(var*mesh.Ne+k)*Np+i] -= integral * quad.wq[q];
				}
			}
		}
	}
}
