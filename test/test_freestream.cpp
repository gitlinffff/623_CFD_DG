#include "surf.hpp"
#include "term2.hpp"
#include "quad.hpp"
#include "problem.hpp"
#include "flux.hpp"
#include "physics.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

namespace {
struct EdgeKey {
    int a;
    int b;
    bool operator<(const EdgeKey& o) const {
        return (a < o.a) || (a == o.a && b < o.b);
    }
};

struct EdgeOwner {
    int elem;
    int face;
    int v0;
    int v1;
};

static double tri_area(double x0, double y0, double x1, double y1, double x2, double y2) {
    return 0.5 * std::abs((x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0));
}

static void outward_normal_and_length(const GriMesh& mesh, int elem, int face, double n_out[2], double& L) {
    int v0 = mesh.E[elem * 3 + (face + 1) % 3];
    int v1 = mesh.E[elem * 3 + (face + 2) % 3];
    double x0 = mesh.V[2 * v0 + 0], y0 = mesh.V[2 * v0 + 1];
    double x1 = mesh.V[2 * v1 + 0], y1 = mesh.V[2 * v1 + 1];
    double ex = x1 - x0;
    double ey = y1 - y0;
    L = std::sqrt(ex * ex + ey * ey);
    if (L < 1e-14) L = 1e-14;
    // For CCW triangle with edge direction v0->v1 along the CCW boundary,
    // outward normal is (ey, -ex)/|e|.
    n_out[0] = ey / L;
    n_out[1] = -ex / L;
}
} // namespace

static void build_rect_mesh_1x1(GriMesh& mesh, int Nx, int Ny) {
    if (Nx < 1) Nx = 1;
    if (Ny < 1) Ny = 1;

    const int npx = Nx + 1;
    const int npy = Ny + 1;
    mesh.Nn = npx * npy;
    mesh.V.resize(mesh.Nn * 2);

    for (int j = 0; j < npy; ++j) {
        for (int i = 0; i < npx; ++i) {
            int vid = j * npx + i;
            mesh.V[2 * vid + 0] = double(i) / double(Nx);
            mesh.V[2 * vid + 1] = double(j) / double(Ny);
        }
    }

    // 2 triangles per quad cell
    mesh.Ne = 2 * Nx * Ny;
    mesh.E.resize(mesh.Ne * 3);
    mesh.Area.resize(mesh.Ne);

    int e = 0;
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            int v00 = j * npx + i;
            int v10 = j * npx + (i + 1);
            int v01 = (j + 1) * npx + i;
            int v11 = (j + 1) * npx + (i + 1);

            // Triangle 1: (v00, v10, v11) CCW
            mesh.E[3 * e + 0] = v00;
            mesh.E[3 * e + 1] = v10;
            mesh.E[3 * e + 2] = v11;
            {
                double x0 = mesh.V[2 * v00 + 0], y0 = mesh.V[2 * v00 + 1];
                double x1 = mesh.V[2 * v10 + 0], y1 = mesh.V[2 * v10 + 1];
                double x2 = mesh.V[2 * v11 + 0], y2 = mesh.V[2 * v11 + 1];
                mesh.Area[e] = tri_area(x0, y0, x1, y1, x2, y2);
            }
            ++e;

            // Triangle 2: (v00, v11, v01) CCW
            mesh.E[3 * e + 0] = v00;
            mesh.E[3 * e + 1] = v11;
            mesh.E[3 * e + 2] = v01;
            {
                double x0 = mesh.V[2 * v00 + 0], y0 = mesh.V[2 * v00 + 1];
                double x1 = mesh.V[2 * v11 + 0], y1 = mesh.V[2 * v11 + 1];
                double x2 = mesh.V[2 * v01 + 0], y2 = mesh.V[2 * v01 + 1];
                mesh.Area[e] = tri_area(x0, y0, x1, y1, x2, y2);
            }
            ++e;
        }
    }

    // Build interior/boundary face tables.
    std::map<EdgeKey, EdgeOwner> edge_first;
    std::vector<int> I2E;
    std::vector<double> In;
    std::vector<double> In_len;

    std::vector<int> B2E;
    std::vector<double> Bn;
    std::vector<double> Bn_len;

    mesh.Bname = {"freestream"};

    for (int elem = 0; elem < mesh.Ne; ++elem) {
        for (int face = 0; face < 3; ++face) {
            int vv0 = mesh.E[elem * 3 + (face + 1) % 3];
            int vv1 = mesh.E[elem * 3 + (face + 2) % 3];
            EdgeKey key{std::min(vv0, vv1), std::max(vv0, vv1)};
            auto it = edge_first.find(key);
            if (it == edge_first.end()) {
                edge_first.emplace(key, EdgeOwner{elem, face, vv0, vv1});
            } else {
                // Interior face between it->second.elem (L) and elem (R)
                int elemL = it->second.elem;
                int faceL = it->second.face;
                int elemR = elem;
                int faceR = face;

                double n_out[2];
                double L;
                outward_normal_and_length(mesh, elemL, faceL, n_out, L);

                I2E.push_back(elemL);
                I2E.push_back(faceL);
                I2E.push_back(elemR);
                I2E.push_back(faceR);
                In.push_back(n_out[0]);
                In.push_back(n_out[1]);
                In_len.push_back(L);

                edge_first.erase(it);
            }
        }
    }

    // Remaining edges are boundary edges
    for (const auto& kv : edge_first) {
        const EdgeOwner& ow = kv.second;
        double n_out[2];
        double L;
        outward_normal_and_length(mesh, ow.elem, ow.face, n_out, L);

        B2E.push_back(ow.elem);
        B2E.push_back(ow.face);
        B2E.push_back(1); // bgroup 1 -> "freestream"
        Bn.push_back(n_out[0]);
        Bn.push_back(n_out[1]);
        Bn_len.push_back(L);
    }

    mesh.num_interior_faces = (int)In_len.size();
    mesh.I2E = std::move(I2E);
    mesh.In = std::move(In);
    mesh.In_len = std::move(In_len);

    mesh.num_boundary_faces = (int)Bn_len.size();
    mesh.B2E = std::move(B2E);
    mesh.Bn = std::move(Bn);
    mesh.Bn_len = std::move(Bn_len);
}

static void write_residual_csv(const GriMesh& mesh, const double* R, int order, const char* path) {
    int Np = (order + 1) * (order + 2) / 2;
    std::ofstream out(path);
    out << "elem,x,y,R_L2\n";

    for (int k = 0; k < mesh.Ne; ++k) {
        int v0 = mesh.E[k * 3 + 0];
        int v1 = mesh.E[k * 3 + 1];
        int v2 = mesh.E[k * 3 + 2];
        double x = (mesh.V[2 * v0 + 0] + mesh.V[2 * v1 + 0] + mesh.V[2 * v2 + 0]) / 3.0;
        double y = (mesh.V[2 * v0 + 1] + mesh.V[2 * v1 + 1] + mesh.V[2 * v2 + 1]) / 3.0;

        double sum = 0.0;
        for (int var = 0; var < 4; ++var) {
            for (int i = 0; i < Np; ++i) {
                double val = R[(var * mesh.Ne + k) * Np + i];
                sum += val * val;
            }
        }
        out << k << "," << x << "," << y << "," << std::sqrt(sum) << "\n";
    }
}

static void write_mesh_and_residual_csv(const GriMesh& mesh, const double* R, int order,
                                        const char* nodes_path, const char* elems_path) {
    int Np = (order + 1) * (order + 2) / 2;

    {
        std::ofstream out(nodes_path);
        out << "vid,x,y\n";
        for (int vid = 0; vid < mesh.Nn; ++vid)
            out << vid << "," << mesh.V[2 * vid + 0] << "," << mesh.V[2 * vid + 1] << "\n";
    }

    {
        std::ofstream out(elems_path);
        out << "elem,v0,v1,v2,cx,cy,R_L2\n";
        for (int k = 0; k < mesh.Ne; ++k) {
            int v0 = mesh.E[k * 3 + 0];
            int v1 = mesh.E[k * 3 + 1];
            int v2 = mesh.E[k * 3 + 2];
            double cx = (mesh.V[2 * v0 + 0] + mesh.V[2 * v1 + 0] + mesh.V[2 * v2 + 0]) / 3.0;
            double cy = (mesh.V[2 * v0 + 1] + mesh.V[2 * v1 + 1] + mesh.V[2 * v2 + 1]) / 3.0;

            double sum = 0.0;
            for (int var = 0; var < 4; ++var) {
                for (int i = 0; i < Np; ++i) {
                    double val = R[(var * mesh.Ne + k) * Np + i];
                    sum += val * val;
                }
            }
            out << k << "," << v0 << "," << v1 << "," << v2 << "," << cx << "," << cy << "," << std::sqrt(sum) << "\n";
        }
    }
}

static void write_solution_csv(const GriMesh& mesh, const double* U, int order,
                               const ProblemParams& params, const char* path) {
    int Np = (order + 1) * (order + 2) / 2;
    std::ofstream out(path);
    out << "elem,cx,cy,rho,u,v,p\n";

    for (int k = 0; k < mesh.Ne; ++k) {
        int v0 = mesh.E[k * 3 + 0];
        int v1 = mesh.E[k * 3 + 1];
        int v2 = mesh.E[k * 3 + 2];
        double cx = (mesh.V[2 * v0 + 0] + mesh.V[2 * v1 + 0] + mesh.V[2 * v2 + 0]) / 3.0;
        double cy = (mesh.V[2 * v0 + 1] + mesh.V[2 * v1 + 1] + mesh.V[2 * v2 + 1]) / 3.0;

        // Freestream init sets all DG DOFs to the same conserved state; use the first DOF.
        double Uk[4];
        for (int var = 0; var < 4; ++var)
            Uk[var] = U[var * mesh.Ne * Np + k * Np + 0];

        double rho, u, v, p, c;
        consToPrim(Uk, params.gammad, rho, u, v, p, c);
        out << k << "," << cx << "," << cy << "," << rho << "," << u << "," << v << "," << p << "\n";
    }
}

int main() {
    const int order = 2;
    const int Np = (order + 1) * (order + 2) / 2;

    GriMesh mesh;
    // 1x1 rectangle domain [0,1]x[0,1], triangulated
    const int Nx = 20;
    const int Ny = 20;
    build_rect_mesh_1x1(mesh, Nx, Ny);

    ProblemParams params;
    const double Mach = 0.2;

    std::vector<double> U(4 * mesh.Ne * Np, 0.0);
    std::vector<double> R(4 * mesh.Ne * Np, 0.0);

    initialize_uniform(U.data(), mesh.Ne, order, params);

    std::vector<BasisEval> phiq = computePhiQ(order);

    // Residual for DG weak form: volume + interior faces + boundary faces.
    addTerm2(mesh, R.data(), order, phiq, U.data(), params);
    addSurfTerm(mesh, R.data(), order, U.data(), params, fluxROE);
    addBndSurfTerm(mesh, R.data(), order, U.data(), params, fluxROE);

    double max_abs = 0.0;
    for (double v : R) max_abs = std::max(max_abs, std::abs(v));

    const double tol = 5e-11; // tight enough for order<=2 on a tiny mesh
    std::cout << "Freestream DG residual max|R| = " << max_abs << "\n";

    write_residual_csv(mesh, R.data(), order, "freestream_residual.csv");
    write_mesh_and_residual_csv(mesh, R.data(), order, "freestream_nodes.csv", "freestream_elems.csv");
    write_solution_csv(mesh, U.data(), order, params, "freestream_solution_elems.csv");
    std::cout << "Wrote: freestream_residual.csv, freestream_nodes.csv, freestream_elems.csv, freestream_solution_elems.csv\n";

    if (max_abs > tol) {
        std::cerr << "FAIL: residual is not near zero (tol=" << tol << ")\n";
        return 1;
    }
    std::cout << "PASS\n";
    return 0;
}

