#include "geometry.hpp"
#include "quad.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>

extern "C" {
    void shapeL(double* xref, int p, double** pphi);
    int gradientL(double* xref, int p, double** pgphi);
}

namespace {

struct ElemGeomView {
    int q;
    int nnode;
    const int* conn;
};

ElemGeomView element_geom_view(const GriMesh& mesh, int elem) {
    if (elem < 0 || elem >= mesh.Ne)
        throw std::out_of_range("element_geom_view: elem out of range");

    if ((int)mesh.ElemGeomOrder.size() == mesh.Ne &&
        (int)mesh.ElemGeomOffset.size() == mesh.Ne + 1) {
        int q = mesh.ElemGeomOrder[elem];
        int nnode = geom_num_nodes(q);
        int off0 = mesh.ElemGeomOffset[elem];
        int off1 = mesh.ElemGeomOffset[elem + 1];
        if (nnode > 0 &&
            off0 >= 0 &&
            off1 >= off0 &&
            off1 <= (int)mesh.ElemGeomConn.size() &&
            off1 - off0 == nnode) {
            return {q, nnode, &mesh.ElemGeomConn[off0]};
        }
    }

    return {1, 3, &mesh.E[elem * 3]};
}

struct AreaCache {
    const GriMesh* mesh;
    std::vector<double> area;
    AreaCache() : mesh(nullptr) {}
};

AreaCache g_area_cache;

void rebuild_area_cache(const GriMesh& mesh) {
    g_area_cache.mesh = &mesh;
    g_area_cache.area.assign(mesh.Ne, 0.0);

    // det(J) for q<=2 geometry is polynomial degree <=2, so order-2 rule is enough.
    QuadratureRule quad = getQuadratureRule(2);

    for (int k = 0; k < mesh.Ne; ++k) {
        double area_k = 0.0;
        for (int q = 0; q < quad.nq; ++q) {
            double xi = quad.xq[2 * q];
            double eta = quad.xq[2 * q + 1];
            GeomEval g;
            eval_geometry(mesh, k, xi, eta, g);
            area_k += quad.wq[q] * std::abs(g.detJ);
        }
        g_area_cache.area[k] = area_k;
    }
}

} // namespace

int geom_num_nodes(int q) {
    if (q < 1) return 0;
    return (q + 1) * (q + 2) / 2;
}

int elem_geom_order(const GriMesh& mesh, int elem) {
    if (elem < 0 || elem >= mesh.Ne) return 1;
    if ((int)mesh.ElemGeomOrder.size() == mesh.Ne)
        return mesh.ElemGeomOrder[elem];
    return 1;
}

void faceRefCoords(int face, double t, double& xi, double& eta) {
    switch (face) {
        case 0: xi = t;       eta = 0.0;     break;
        case 1: xi = 1.0 - t; eta = t;       break;
        case 2: xi = 0.0;     eta = 1.0 - t; break;
        default: xi = t;      eta = 0.0;     break;
    }
}

void faceRefDerivatives(int face, double& dxi_dt, double& deta_dt) {
    switch (face) {
        case 0: dxi_dt = 1.0;  deta_dt = 0.0;  break;
        case 1: dxi_dt = -1.0; deta_dt = 1.0;  break;
        case 2: dxi_dt = 0.0;  deta_dt = -1.0; break;
        default: dxi_dt = 1.0; deta_dt = 0.0;  break;
    }
}

void eval_geometry(const GriMesh& mesh, int elem, double xi, double eta, GeomEval& g) {
    ElemGeomView eg = element_geom_view(mesh, elem);

    if (eg.q == 1) {
        int n0 = eg.conn[0];
        int n1 = eg.conn[1];
        int n2 = eg.conn[2];
        double x0 = mesh.V[n0 * 2 + 0], y0 = mesh.V[n0 * 2 + 1];
        double x1 = mesh.V[n1 * 2 + 0], y1 = mesh.V[n1 * 2 + 1];
        double x2 = mesh.V[n2 * 2 + 0], y2 = mesh.V[n2 * 2 + 1];

        g.x = x0 + xi * (x1 - x0) + eta * (x2 - x0);
        g.y = y0 + xi * (y1 - y0) + eta * (y2 - y0);
        g.J(0, 0) = x1 - x0; g.J(0, 1) = x2 - x0;
        g.J(1, 0) = y1 - y0; g.J(1, 1) = y2 - y0;
        g.detJ = g.J.determinant();
        return;
    }

    if (eg.q < 1 || eg.nnode <= 0)
        throw std::runtime_error("eval_geometry: unsupported geometric order");

    thread_local double* phi = nullptr;
    thread_local double* gphi = nullptr;
    double xref[2] = {xi, eta};
    shapeL(xref, eg.q, &phi);
    if (gradientL(xref, eg.q, &gphi) != 0)
        throw std::runtime_error("eval_geometry: gradientL failed");

    g.x = 0.0;
    g.y = 0.0;
    g.J.setZero();
    for (int i = 0; i < eg.nnode; ++i) {
        int nid = eg.conn[i];
        double x = mesh.V[nid * 2 + 0];
        double y = mesh.V[nid * 2 + 1];
        double dphi_dxi = gphi[i];
        double dphi_deta = gphi[i + eg.nnode];

        g.x += x * phi[i];
        g.y += y * phi[i];
        g.J(0, 0) += x * dphi_dxi;
        g.J(0, 1) += x * dphi_deta;
        g.J(1, 0) += y * dphi_dxi;
        g.J(1, 1) += y * dphi_deta;
    }
    g.detJ = g.J.determinant();
}

bool phys_to_ref(const GriMesh& mesh, int elem, double x, double y, double& xi, double& eta) {
    int q = elem_geom_order(mesh, elem);

    // Start from linear-corner mapping.
    int n0 = mesh.E[elem * 3 + 0];
    int n1 = mesh.E[elem * 3 + 1];
    int n2 = mesh.E[elem * 3 + 2];
    Eigen::Matrix2d J0;
    J0(0, 0) = mesh.V[n1 * 2 + 0] - mesh.V[n0 * 2 + 0];
    J0(0, 1) = mesh.V[n2 * 2 + 0] - mesh.V[n0 * 2 + 0];
    J0(1, 0) = mesh.V[n1 * 2 + 1] - mesh.V[n0 * 2 + 1];
    J0(1, 1) = mesh.V[n2 * 2 + 1] - mesh.V[n0 * 2 + 1];
    double det0 = J0.determinant();
    if (std::abs(det0) < 1e-14) return false;
    Eigen::Vector2d d0;
    d0 << x - mesh.V[n0 * 2 + 0], y - mesh.V[n0 * 2 + 1];
    Eigen::Vector2d ref0 = J0.inverse() * d0;
    xi = ref0(0);
    eta = ref0(1);

    if (q == 1) return true;

    const int max_iter = 30;
    const double tol = 1e-12;
    for (int it = 0; it < max_iter; ++it) {
        GeomEval g;
        eval_geometry(mesh, elem, xi, eta, g);
        double rx = x - g.x;
        double ry = y - g.y;
        double rnorm = std::sqrt(rx * rx + ry * ry);
        if (rnorm < tol) return true;

        if (std::abs(g.detJ) < 1e-14) return false;
        Eigen::Vector2d rhs;
        rhs << rx, ry;
        Eigen::Vector2d delta = g.J.inverse() * rhs;
        xi += delta(0);
        eta += delta(1);

        if (!std::isfinite(xi) || !std::isfinite(eta)) return false;
    }

    GeomEval g;
    eval_geometry(mesh, elem, xi, eta, g);
    double rx = x - g.x;
    double ry = y - g.y;
    double rnorm = std::sqrt(rx * rx + ry * ry);
    return (rnorm < 1e-8);
}

double face_metric_normal_point(const GriMesh& mesh, int elem, int face, double t,
                                double n_out[2], double& x, double& y) {
    double xi, eta;
    faceRefCoords(face, t, xi, eta);
    GeomEval g;
    eval_geometry(mesh, elem, xi, eta, g);

    double dxi_dt, deta_dt;
    faceRefDerivatives(face, dxi_dt, deta_dt);

    // Tangent = dx/dt = J * d(xi,eta)/dt
    double tx = g.J(0, 0) * dxi_dt + g.J(0, 1) * deta_dt;
    double ty = g.J(1, 0) * dxi_dt + g.J(1, 1) * deta_dt;
    double ds = std::sqrt(tx * tx + ty * ty);
    if (ds < 1e-14) ds = 1e-14;

    // For CCW element and boundary parametrised along local face ordering, (ty,-tx) is outward.
    n_out[0] = ty / ds;
    n_out[1] = -tx / ds;

    x = g.x;
    y = g.y;
    return ds;
}

const std::vector<double>& get_geometric_element_areas(const GriMesh& mesh) {
    if (g_area_cache.mesh != &mesh || (int)g_area_cache.area.size() != mesh.Ne)
        rebuild_area_cache(mesh);
    return g_area_cache.area;
}

double get_geometric_element_area(const GriMesh& mesh, int elem) {
    const std::vector<double>& area = get_geometric_element_areas(mesh);
    if (elem < 0 || elem >= (int)area.size())
        throw std::out_of_range("get_geometric_element_area: elem out of range");
    return area[elem];
}
