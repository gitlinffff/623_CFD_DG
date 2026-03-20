#include "mm.hpp"
#include "geometry.hpp"
#include <cmath>
#include <stdexcept>

extern "C" {
    void shapeL(double* xref, int p, double** pphi);
}

namespace {

Eigen::MatrixXd assemble_element_mass_matrix(const GriMesh& mesh, int elem, int order,
                                             const QuadratureRule& quad,
                                             const std::vector<BasisEval>& phiq) {
    const int Np = (order + 1) * (order + 2) / 2;
    Eigen::MatrixXd Mk = Eigen::MatrixXd::Zero(Np, Np);

    const bool affine_geom = (elem_geom_order(mesh, elem) == 1);
    double detJ_affine = 0.0;
    if (affine_geom) {
        GeomEval g0;
        eval_geometry(mesh, elem, 1.0 / 3.0, 1.0 / 3.0, g0);
        detJ_affine = std::abs(g0.detJ);
    }

    for (int q = 0; q < quad.nq; ++q) {
        double detJ = 0.0;
        if (affine_geom) {
            detJ = detJ_affine;
        } else {
            double xi = quad.xq[2 * q];
            double eta = quad.xq[2 * q + 1];
            GeomEval gq;
            eval_geometry(mesh, elem, xi, eta, gq);
            detJ = std::abs(gq.detJ);
        }
        if (detJ < 1e-30) detJ = 1e-30;

        for (int i = 0; i < Np; ++i) {
            double phi_i = phiq[i * quad.nq + q].phi;
            for (int j = 0; j < Np; ++j) {
                double phi_j = phiq[j * quad.nq + q].phi;
                Mk(i, j) += quad.wq[q] * detJ * phi_i * phi_j;
            }
        }
    }
    return Mk;
}

struct MassInvCache {
    const GriMesh* mesh;
    int order;
    std::vector<Eigen::MatrixXd> Minv_elem;
    MassInvCache() : mesh(nullptr), order(-1) {}
};

MassInvCache g_mass_inv_cache;

void rebuild_mass_inverse_cache(const GriMesh& mesh, int order) {
    const int Np = (order + 1) * (order + 2) / 2;
    QuadratureRule quad = getQuadratureRule(order);
    std::vector<BasisEval> phiq = computePhiQ(order);

    g_mass_inv_cache.mesh = &mesh;
    g_mass_inv_cache.order = order;
    g_mass_inv_cache.Minv_elem.assign(mesh.Ne, Eigen::MatrixXd::Zero(Np, Np));

    for (int k = 0; k < mesh.Ne; ++k) {
        Eigen::MatrixXd Mk = assemble_element_mass_matrix(mesh, k, order, quad, phiq);
        g_mass_inv_cache.Minv_elem[k] = Mk.inverse();
    }
}

} // namespace

Eigen::Matrix2d Jacobian(const GriMesh& mesh, int elem) {
    // Legacy/utility corner Jacobian for the affine triangle built from corner nodes.
    Eigen::Matrix2d J;
    int v1_idx = mesh.E[elem * 3 + 0];
    int v2_idx = mesh.E[elem * 3 + 1];
    int v3_idx = mesh.E[elem * 3 + 2];

    double x1 = mesh.V[v1_idx * 2 + 0], y1 = mesh.V[v1_idx * 2 + 1];
    double x2 = mesh.V[v2_idx * 2 + 0], y2 = mesh.V[v2_idx * 2 + 1];
    double x3 = mesh.V[v3_idx * 2 + 0], y3 = mesh.V[v3_idx * 2 + 1];

    J << x2 - x1, x3 - x1,
         y2 - y1, y3 - y1;
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
            for (int j = 0; j < Np; ++j)
                M_ref(i, j) += quad.wq[q] * phi[i] * phi[j];
    }

    free(phi);
    return M_ref;
}

Eigen::SparseMatrix<double> computeGlobalMassMatrix(const GriMesh& mesh, int order) {
    int Np = (order + 1) * (order + 2) / 2;
    int totalSize = mesh.Ne * Np;

    Eigen::SparseMatrix<double> M(totalSize, totalSize);
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve((size_t)mesh.Ne * (size_t)Np * (size_t)Np);

    QuadratureRule quad = getQuadratureRule(order);
    std::vector<BasisEval> phiq = computePhiQ(order);

    for (int k = 0; k < mesh.Ne; ++k) {
        Eigen::MatrixXd Mk = assemble_element_mass_matrix(mesh, k, order, quad, phiq);
        int blockStart = k * Np;
        for (int i = 0; i < Np; ++i) {
            for (int j = 0; j < Np; ++j) {
                tripletList.push_back(Eigen::Triplet<double>(blockStart + i, blockStart + j, Mk(i, j)));
            }
        }
    }

    M.setFromTriplets(tripletList.begin(), tripletList.end());
    M.makeCompressed();
    return M;
}

void applyInverseMassMatrix(const GriMesh& mesh, double* R, int order) {
    int Np = (order + 1) * (order + 2) / 2;

    if (g_mass_inv_cache.mesh != &mesh || g_mass_inv_cache.order != order ||
        (int)g_mass_inv_cache.Minv_elem.size() != mesh.Ne) {
        rebuild_mass_inverse_cache(mesh, order);
    }

    const std::vector<Eigen::MatrixXd>& Minv_elem = g_mass_inv_cache.Minv_elem;

    #pragma omp parallel for schedule(static)
    for (int k = 0; k < mesh.Ne; ++k) {
        Eigen::VectorXd R_elem(Np);
        for (int var = 0; var < 4; ++var) {
            int base = (var * mesh.Ne + k) * Np;
            for (int i = 0; i < Np; ++i)
                R_elem(i) = R[base + i];
            Eigen::VectorXd R_new = Minv_elem[k] * R_elem;
            for (int i = 0; i < Np; ++i)
                R[base + i] = R_new(i);
        }
    }
}
