#include "write_vtu.hpp"
#include "physics.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sys/stat.h>

extern "C" {
    void shapeL(double* xref, int p, double** pphi);
}

static void make_dirs(const std::string& path)
{
    for (size_t i = 1; i <= path.size(); ++i) {
        if (i == path.size() || path[i] == '/') {
            std::string sub = path.substr(0, i);
            mkdir(sub.c_str(), 0755);
        }
    }
}

static void quadratic_edge_point(const GriMesh& mesh, int n0, int n1, int nm, double t,
                                 double& x, double& y)
{
    const double x0 = mesh.V[n0 * 2 + 0], y0 = mesh.V[n0 * 2 + 1];
    const double x1 = mesh.V[n1 * 2 + 0], y1 = mesh.V[n1 * 2 + 1];
    const double xm = mesh.V[nm * 2 + 0], ym = mesh.V[nm * 2 + 1];

    const double N0 = (1.0 - t) * (1.0 - 2.0 * t);
    const double Nm = 4.0 * t * (1.0 - t);
    const double N1 = t * (2.0 * t - 1.0);

    x = N0 * x0 + Nm * xm + N1 * x1;
    y = N0 * y0 + Nm * ym + N1 * y1;
}

static void linear_edge_point(const GriMesh& mesh, int n0, int n1, double t,
                              double& x, double& y)
{
    const double x0 = mesh.V[n0 * 2 + 0], y0 = mesh.V[n0 * 2 + 1];
    const double x1 = mesh.V[n1 * 2 + 0], y1 = mesh.V[n1 * 2 + 1];
    x = (1.0 - t) * x0 + t * x1;
    y = (1.0 - t) * y0 + t * y1;
}

static int visualization_subdivisions(int order)
{
    switch (order) {
        case 0: return 3;
        case 1: return 4;
        case 2: return 6;
        case 3: return 8;
        default:
            return std::max(order + 1, 8);
    }
}

static std::string lower_copy(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

static bool has_token(const std::string& s, const char* token)
{
    return s.find(token) != std::string::npos;
}

static bool is_top_wall_name(const std::string& bname_raw)
{
    std::string bname = lower_copy(bname_raw);
    if (bname == "bgroup2") return true;
    return has_token(bname, "top") || has_token(bname, "upper");
}

static bool is_bottom_wall_name(const std::string& bname_raw)
{
    std::string bname = lower_copy(bname_raw);
    if (bname == "bgroup6") return true;
    return has_token(bname, "bottom") || has_token(bname, "lower");
}

void write_solution_vtu(const GriMesh& mesh, const double* U, int order,
                        const ProblemParams& params, const std::string& filepath)
{
    const int Np     = (order + 1) * (order + 2) / 2;
    const int n_ref  = visualization_subdivisions(order);
    const double h   = 1.0 / n_ref;

    std::vector<std::vector<int>> vtx_flat(n_ref + 1,
                                           std::vector<int>(n_ref + 1, -1));
    std::vector<double> sv_xi, sv_eta;

    for (int j = 0; j <= n_ref; ++j)
        for (int i = 0; i <= n_ref - j; ++i) {
            vtx_flat[i][j] = (int)sv_xi.size();
            sv_xi.push_back(i * h);
            sv_eta.push_back(j * h);
        }

    int n_sv = (int)sv_xi.size();

    std::vector<int> tv0, tv1, tv2;
    for (int j = 0; j < n_ref; ++j) {
        for (int i = 0; i < n_ref - j; ++i) {

            tv0.push_back(vtx_flat[i    ][j    ]);
            tv1.push_back(vtx_flat[i + 1][j    ]);
            tv2.push_back(vtx_flat[i    ][j + 1]);

            if (i + j + 1 < n_ref) {
                tv0.push_back(vtx_flat[i + 1][j    ]);
                tv1.push_back(vtx_flat[i + 1][j + 1]);
                tv2.push_back(vtx_flat[i    ][j + 1]);
            }
        }
    }
    int n_st = (int)tv0.size();

    long total_pts   = (long)mesh.Ne * n_sv;
    long total_cells = (long)mesh.Ne * n_st;

    std::vector<int> curved_mid(mesh.Ne * 3, -1);
    std::vector<int> wall_face_mark(mesh.Ne * 3, 0);
    if ((int)mesh.BedgeNodeOffset.size() == mesh.num_boundary_faces + 1) {
        for (int ib = 0; ib < mesh.num_boundary_faces; ++ib) {
            int start = mesh.BedgeNodeOffset[(size_t)ib];
            int nnode = mesh.BedgeNodeOffset[(size_t)ib + 1] - start;
            int elem = mesh.B2E[3 * ib + 0];
            int face = mesh.B2E[3 * ib + 1];
            int bgroup = mesh.B2E[3 * ib + 2];
            const int local_edge_1based = face + 1;
            const int* tri = &mesh.E[elem * 3];
            int n0 = tri[face];
            int n1 = tri[(face + 1) % 3];

            if (bgroup >= 1 && (size_t)bgroup <= mesh.Bname.size()) {
                const std::string& bname = mesh.Bname[(size_t)bgroup - 1];
                if (is_top_wall_name(bname)) {
                    wall_face_mark[elem * 3 + face] = local_edge_1based;
                } else if (is_bottom_wall_name(bname)) {
                    wall_face_mark[elem * 3 + face] = -local_edge_1based;
                }
            }

            if (nnode != 3) continue;

            if (mesh.BedgeNodes[(size_t)start + 0] == n0 &&
                mesh.BedgeNodes[(size_t)start + 1] == n1) {
                curved_mid[elem * 3 + face] = mesh.BedgeNodes[(size_t)start + 2];
            }
        }
    }

    std::vector<std::vector<double>> phi_sv(n_sv, std::vector<double>(Np));
    {
        double* buf = nullptr;
        for (int v = 0; v < n_sv; ++v) {
            double xref[2] = { sv_xi[v], sv_eta[v] };
            shapeL(xref, order, &buf);
            for (int j = 0; j < Np; ++j)
                phi_sv[v][j] = buf[j];
        }
        if (buf) free(buf);
    }

    std::vector<double> f_rho  (mesh.Ne * n_sv);
    std::vector<double> f_u    (mesh.Ne * n_sv);
    std::vector<double> f_v    (mesh.Ne * n_sv);
    std::vector<double> f_p    (mesh.Ne * n_sv);
    std::vector<double> f_mach (mesh.Ne * n_sv);
    std::vector<double> f_entr (mesh.Ne * n_sv);
    std::vector<double> f_wall_marker(mesh.Ne * n_sv, 0.0);

    const double s_ref = params.a0 * params.a0 / (params.gammad *
                         std::pow(params.rho0, params.gammad - 1.0));

    for (int k = 0; k < mesh.Ne; ++k) {
        for (int v = 0; v < n_sv; ++v) {
            const double xi  = sv_xi[v];
            const double eta = sv_eta[v];
            const double L0 = 1.0 - xi - eta;
            const double L1 = xi;
            const double L2 = eta;
            const double tol_edge = 1e-12;

            int wall_marker = 0;
            if (wall_face_mark[k * 3 + 0] != 0 && std::fabs(L2) <= tol_edge) {
                wall_marker = wall_face_mark[k * 3 + 0];
            }
            if (wall_face_mark[k * 3 + 1] != 0 && std::fabs(L0) <= tol_edge) {
                if (wall_marker == 0) wall_marker = wall_face_mark[k * 3 + 1];
            }
            if (wall_face_mark[k * 3 + 2] != 0 && std::fabs(L1) <= tol_edge) {
                if (wall_marker == 0) wall_marker = wall_face_mark[k * 3 + 2];
            }

            double Uk[4] = {0, 0, 0, 0};
            for (int var = 0; var < 4; ++var)
                for (int j = 0; j < Np; ++j)
                    Uk[var] += U[var * mesh.Ne * Np + k * Np + j] * phi_sv[v][j];

            double rho, uu, vv, p, c;
            consToPrim(Uk, params.gammad, rho, uu, vv, p, c);

            double spd  = std::sqrt(uu * uu + vv * vv);
            double mach = (c > 0) ? spd / c : 0.0;
            double entr = (rho > 0) ? (p / std::pow(rho, params.gammad)) / s_ref
                                     : 1.0;

            int idx = k * n_sv + v;
            f_rho [idx] = rho;
            f_u   [idx] = uu;
            f_v   [idx] = vv;
            f_p   [idx] = p;
            f_mach[idx] = mach;
            f_entr[idx] = entr;
            f_wall_marker[idx] = (double)wall_marker;
        }
    }

    {
        size_t slash = filepath.find_last_of('/');
        if (slash != std::string::npos)
            make_dirs(filepath.substr(0, slash));
    }

    std::ofstream out(filepath.c_str());
    if (!out) {
        std::cerr << "[write_vtu] Cannot open file: " << filepath << "\n";
        return;
    }

    out << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\""
        <<         " byte_order=\"LittleEndian\">\n"
        << "  <UnstructuredGrid>\n"
        << "    <Piece NumberOfPoints=\"" << total_pts
        << "\" NumberOfCells=\"" << total_cells << "\">\n";

    out << "      <Points>\n"
        << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\""
        <<                   " format=\"ascii\">\n";
    for (int k = 0; k < mesh.Ne; ++k) {
        int n0 = mesh.E[k * 3 + 0];
        int n1 = mesh.E[k * 3 + 1];
        int n2 = mesh.E[k * 3 + 2];
        double x0 = mesh.V[n0 * 2], y0 = mesh.V[n0 * 2 + 1];
        double x1 = mesh.V[n1 * 2], y1 = mesh.V[n1 * 2 + 1];
        double x2 = mesh.V[n2 * 2], y2 = mesh.V[n2 * 2 + 1];

        for (int v = 0; v < n_sv; ++v) {
            double xi  = sv_xi[v], eta = sv_eta[v];
            double xp  = x0 + xi * (x1 - x0) + eta * (x2 - x0);
            double yp  = y0 + xi * (y1 - y0) + eta * (y2 - y0);

            const double L0 = 1.0 - xi - eta;
            const double L1 = xi;
            const double L2 = eta;
            const double eps = 1e-14;

            auto add_face_correction = [&](int face, int mid_node) {
                if (mid_node < 0) return;

                int a, b;
                double den, t;
                if (face == 0) {
                    a = n0; b = n1;
                    den = L0 + L1;
                    if (den <= eps) return;
                    t = L1 / den;
                } else if (face == 1) {
                    a = n1; b = n2;
                    den = L1 + L2;
                    if (den <= eps) return;
                    t = L2 / den;
                } else {
                    a = n2; b = n0;
                    den = L2 + L0;
                    if (den <= eps) return;
                    t = L0 / den;
                }

                t = std::min(1.0, std::max(0.0, t));
                double xc, yc, xl, yl;
                quadratic_edge_point(mesh, a, b, mid_node, t, xc, yc);
                linear_edge_point(mesh, a, b, t, xl, yl);
                xp += den * (xc - xl);
                yp += den * (yc - yl);
            };

            add_face_correction(0, curved_mid[k * 3 + 0]);
            add_face_correction(1, curved_mid[k * 3 + 1]);
            add_face_correction(2, curved_mid[k * 3 + 2]);

            out << "          " << xp << " " << yp << " 0\n";
        }
    }
    out << "        </DataArray>\n"
        << "      </Points>\n";

    out << "      <Cells>\n";

    out << "        <DataArray type=\"Int64\" Name=\"connectivity\""
        <<                   " format=\"ascii\">\n";
    for (int k = 0; k < mesh.Ne; ++k) {
        long base = (long)k * n_sv;
        for (int t = 0; t < n_st; ++t)
            out << "          " << (base + tv0[t])
                << " "          << (base + tv1[t])
                << " "          << (base + tv2[t]) << "\n";
    }
    out << "        </DataArray>\n";

    out << "        <DataArray type=\"Int64\" Name=\"offsets\""
        <<                   " format=\"ascii\">\n";
    for (long c = 1; c <= total_cells; ++c)
        out << "          " << c * 3 << "\n";
    out << "        </DataArray>\n";

    out << "        <DataArray type=\"UInt8\" Name=\"types\""
        <<                   " format=\"ascii\">\n";
    for (long c = 0; c < total_cells; ++c)
        out << "          5\n";
    out << "        </DataArray>\n";

    out << "      </Cells>\n";

    out << "      <PointData>\n";

    auto write_scalar = [&](const char* name, const std::vector<double>& field) {
        out << "        <DataArray type=\"Float64\" Name=\"" << name
            << "\" format=\"ascii\">\n";
        for (int k = 0; k < mesh.Ne; ++k)
            for (int v = 0; v < n_sv; ++v)
                out << "          " << field[k * n_sv + v] << "\n";
        out << "        </DataArray>\n";
    };

    write_scalar("rho",     f_rho);
    write_scalar("u",       f_u);
    write_scalar("v",       f_v);
    write_scalar("p",       f_p);
    write_scalar("Mach",    f_mach);
    write_scalar("entropy", f_entr);
    write_scalar("wall_marker", f_wall_marker);

    out << "        <DataArray type=\"Float64\" Name=\"velocity\""
        <<                   " NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int k = 0; k < mesh.Ne; ++k)
        for (int v = 0; v < n_sv; ++v)
            out << "          " << f_u[k * n_sv + v]
                << " "          << f_v[k * n_sv + v]
                << " 0\n";
    out << "        </DataArray>\n";

    out << "      </PointData>\n";

    out << "    </Piece>\n"
        << "  </UnstructuredGrid>\n"
        << "</VTKFile>\n";

    out.close();
    std::cout << "[write_vtu] Solution written to: " << filepath << "\n";
}
