#include "write_vtu.hpp"
#include "physics.hpp"
#include "geometry.hpp"
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

/* ---------- helpers --------------------------------------------------- */

static void make_dirs(const std::string& path)
{
    for (size_t i = 1; i <= path.size(); ++i) {
        if (i == path.size() || path[i] == '/') {
            std::string sub = path.substr(0, i);
            mkdir(sub.c_str(), 0755);   /* ignore EEXIST */
        }
    }
}

enum class WallKind {
    NONE,
    UPPER,
    LOWER
};

static std::string to_lower(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

static bool has_token(const std::string& s, const char* token)
{
    return s.find(token) != std::string::npos;
}

static WallKind classify_wall_group(const std::string& name_raw)
{
    std::string name = to_lower(name_raw);

    // Project-specific blade groups.
    if (name == "bgroup6") return WallKind::UPPER;
    if (name == "bgroup2") return WallKind::LOWER;

    // Generic naming fallback.
    if (has_token(name, "upper") && has_token(name, "wall")) return WallKind::UPPER;
    if (has_token(name, "lower") && has_token(name, "wall")) return WallKind::LOWER;

    return WallKind::NONE;
}

/* ---------- main writer ----------------------------------------------- */

void write_solution_vtu(const GriMesh& mesh, const double* U, int order,
                        const ProblemParams& params, const std::string& filepath,
                        bool log_message)
{
    const int Np     = (order + 1) * (order + 2) / 2;
    const int n_ref  = std::max(1, 2 * order);
    const double h   = 1.0 / n_ref;

    /* ---- Build uniform sub-grid on reference triangle ----------------
     *  Vertices at (i*h, j*h) for i,j >= 0, i+j <= n_ref.
     *  Index stored in vtx_flat[i][j].
     */
    std::vector<std::vector<int>> vtx_flat(n_ref + 1,
                                           std::vector<int>(n_ref + 1, -1));
    std::vector<double> sv_xi, sv_eta;

    for (int j = 0; j <= n_ref; ++j)
        for (int i = 0; i <= n_ref - j; ++i) {
            vtx_flat[i][j] = (int)sv_xi.size();
            sv_xi.push_back(i * h);
            sv_eta.push_back(j * h);
        }

    int n_sv = (int)sv_xi.size();   /* (n_ref+1)*(n_ref+2)/2 */

    /* ---- Build sub-triangle connectivity (CCW orientation) ----------- */
    std::vector<int> tv0, tv1, tv2;
    for (int j = 0; j < n_ref; ++j) {
        for (int i = 0; i < n_ref - j; ++i) {
            /* upward-pointing */
            tv0.push_back(vtx_flat[i    ][j    ]);
            tv1.push_back(vtx_flat[i + 1][j    ]);
            tv2.push_back(vtx_flat[i    ][j + 1]);
            /* downward-pointing (exists when i+j+1 < n_ref) */
            if (i + j + 1 < n_ref) {
                tv0.push_back(vtx_flat[i + 1][j    ]);
                tv1.push_back(vtx_flat[i + 1][j + 1]);
                tv2.push_back(vtx_flat[i    ][j + 1]);
            }
        }
    }
    int n_st = (int)tv0.size();   /* n_ref^2 */

    long total_pts   = (long)mesh.Ne * n_sv;
    long total_cells = (long)mesh.Ne * n_st;

    // Element-local wall marker per face:
    //  0 = non-wall, +{1,2,3}=upper face index, -{1,2,3}=lower face index.
    std::vector<int> elem_wall_face_marker(mesh.Ne * 3, 0);
    for (int ib = 0; ib < mesh.num_boundary_faces; ++ib) {
        int elem = mesh.B2E[3 * ib + 0];
        int face = mesh.B2E[3 * ib + 1];
        int bgrp = mesh.B2E[3 * ib + 2];

        std::string bname = "unknown";
        if (bgrp >= 1 && (size_t)bgrp <= mesh.Bname.size())
            bname = mesh.Bname[(size_t)bgrp - 1];

        WallKind kind = classify_wall_group(bname);
        if (kind == WallKind::NONE) continue;
        int marker = (kind == WallKind::UPPER ? +1 : -1) * (face + 1);
        elem_wall_face_marker[elem * 3 + face] = marker;
    }

    /* ---- Precompute Lagrange basis at sub-vertices (same for all elements) */
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

    /* ---- Evaluate primitive fields at every (element, sub-vertex) ---- */
    /* Layout: field[k * n_sv + v] */
    std::vector<double> f_rho  (mesh.Ne * n_sv);
    std::vector<double> f_u    (mesh.Ne * n_sv);
    std::vector<double> f_v    (mesh.Ne * n_sv);
    std::vector<double> f_p    (mesh.Ne * n_sv);
    std::vector<double> f_mach (mesh.Ne * n_sv);
    std::vector<double> f_entr (mesh.Ne * n_sv);   /* s/s_ref = (p/rho^g) / s_ref */
    std::vector<int>    f_wall_marker(mesh.Ne * n_sv, 0);

    /* Reference entropy from freestream */
    const double s_ref = params.a0 * params.a0 / (params.gammad *
                         std::pow(params.rho0, params.gammad - 1.0));

    for (int k = 0; k < mesh.Ne; ++k) {
        for (int v = 0; v < n_sv; ++v) {
            /* Interpolate conservative state at this sub-vertex */
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

            // Mark wall points using reference-face location.
            const double xi = sv_xi[v];
            const double eta = sv_eta[v];
            const double tol = 1e-12;
            int marker = 0;
            for (int face = 0; face < 3; ++face) {
                int face_marker = elem_wall_face_marker[k * 3 + face];
                if (face_marker == 0) continue;

                bool on_face = false;
                if (face == 0) on_face = std::abs(eta) <= tol;
                if (face == 1) on_face = std::abs(xi + eta - 1.0) <= tol;
                if (face == 2) on_face = std::abs(xi) <= tol;
                if (on_face) {
                    marker = face_marker;
                    break;
                }
            }
            f_wall_marker[idx] = marker;
        }
    }

    /* ---- Create output directory and open file ----------------------- */
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

    /* ---- VTK XML header --------------------------------------------- */
    out << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\""
        <<         " byte_order=\"LittleEndian\">\n"
        << "  <UnstructuredGrid>\n"
        << "    <Piece NumberOfPoints=\"" << total_pts
        << "\" NumberOfCells=\"" << total_cells << "\">\n";

    /* ---- Points (physical coordinates) ------------------------------ */
    out << "      <Points>\n"
        << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\""
        <<                   " format=\"ascii\">\n";
    for (int k = 0; k < mesh.Ne; ++k) {
        for (int v = 0; v < n_sv; ++v) {
            double xi  = sv_xi[v], eta = sv_eta[v];
            GeomEval g;
            eval_geometry(mesh, k, xi, eta, g);
            double xp  = g.x;
            double yp  = g.y;
            out << "          " << xp << " " << yp << " 0\n";
        }
    }
    out << "        </DataArray>\n"
        << "      </Points>\n";

    /* ---- Cells ------------------------------------------------------ */
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
        out << "          5\n";   /* VTK_TRIANGLE */
    out << "        </DataArray>\n";

    out << "      </Cells>\n";

    /* ---- PointData fields ------------------------------------------- */
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

    out << "        <DataArray type=\"Int32\" Name=\"wall_marker\" format=\"ascii\">\n";
    for (int k = 0; k < mesh.Ne; ++k)
        for (int v = 0; v < n_sv; ++v)
            out << "          " << f_wall_marker[k * n_sv + v] << "\n";
    out << "        </DataArray>\n";

    /* Velocity as a vector (3-component: u, v, 0) */
    out << "        <DataArray type=\"Float64\" Name=\"velocity\""
        <<                   " NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int k = 0; k < mesh.Ne; ++k)
        for (int v = 0; v < n_sv; ++v)
            out << "          " << f_u[k * n_sv + v]
                << " "          << f_v[k * n_sv + v]
                << " 0\n";
    out << "        </DataArray>\n";

    out << "      </PointData>\n";

    /* ---- Footer ----------------------------------------------------- */
    out << "    </Piece>\n"
        << "  </UnstructuredGrid>\n"
        << "</VTKFile>\n";

    out.close();
    if (log_message)
        std::cout << "[write_vtu] Solution written to: " << filepath << "\n";
}
