#include "term2.hpp"
#include "physics.hpp"
#include "geometry.hpp"
#include <cmath>

void addTerm2(const GriMesh& mesh, double* R, int order,
              const std::vector<BasisEval>& phiq, const double* U,
              const ProblemParams& params) {
    int Np = (order + 1) * (order + 2) / 2;
    QuadratureRule quad = getQuadratureRule(order);

    // Each element k writes exclusively to R[(var*Ne+k)*Np+i], so no race conditions.
    #pragma omp parallel for schedule(static)
    for (int k = 0; k < mesh.Ne; ++k) {
        const bool affine_geom = (elem_geom_order(mesh, k) == 1);
        Eigen::Matrix2d J_inv_affine;
        double detJ_affine = 0.0;
        if (affine_geom) {
            GeomEval g0;
            eval_geometry(mesh, k, 1.0 / 3.0, 1.0 / 3.0, g0);
            detJ_affine = std::abs(g0.detJ);
            if (detJ_affine < 1e-30) detJ_affine = 1e-30;
            J_inv_affine = g0.J.inverse();
        }

        for (int q = 0; q < quad.nq; ++q) {
            Eigen::Matrix2d J_inv;
            double detJ = 0.0;
            if (affine_geom) {
                J_inv = J_inv_affine;
                detJ = detJ_affine;
            } else {
                double xi = quad.xq[2 * q];
                double eta = quad.xq[2 * q + 1];
                GeomEval gq;
                eval_geometry(mesh, k, xi, eta, gq);
                detJ = std::abs(gq.detJ);
                if (detJ < 1e-30) detJ = 1e-30;
                J_inv = gq.J.inverse();
            }

            // Interpolate conservative state at this quad point.
            Eigen::Vector4d uq = Eigen::Vector4d::Zero();
            for (int i = 0; i < Np; ++i) {
                double phi_val = phiq[i * quad.nq + q].phi;
                for (int var = 0; var < 4; ++var) {
                    uq(var) += U[var * mesh.Ne * Np + k * Np + i] * phi_val;
                }
            }

            // Physical flux at this quadrature point.
            Eigen::Matrix<double, 4, 2> phyF = vec_phyFlux(uq.data(), params.gammad);

            // Accumulate weak-form volume contribution.
            for (int i = 0; i < Np; ++i) {
                Eigen::RowVector2d gradPhi_ref;
                gradPhi_ref << phiq[i * quad.nq + q].dphi_dxi,
                               phiq[i * quad.nq + q].dphi_deta;
                Eigen::RowVector2d gradPhi_phys = gradPhi_ref * J_inv;

                for (int var = 0; var < 4; ++var) {
                    double integral = gradPhi_phys(0) * phyF(var, 0) +
                                      gradPhi_phys(1) * phyF(var, 1);
                    R[(var * mesh.Ne + k) * Np + i] -= integral * detJ * quad.wq[q];
                }
            }
        }
    }
}
