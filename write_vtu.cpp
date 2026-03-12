#include "write_vtu.hpp"
#include "physics.hpp"
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

/* ---------- main writer ----------------------------------------------- */

void write_solution_vtu(const GriMesh& mesh, const double* U, int order,
                        const ProblemParams& params, const std::string& filepath)
{
    const int Np     = (order + 1) * (order + 2) / 2;
    const int n_ref  = order + 1;           /* subdivision level */
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
        int n0 = mesh.E[k * 3 + 0];
        int n1 = mesh.E[k * 3 + 1];
        int n2 = mesh.E[k * 3 + 2];
        double x0 = mesh.V[n0 * 2], y0 = mesh.V[n0 * 2 + 1];
        double x1 = mesh.V[n1 * 2], y1 = mesh.V[n1 * 2 + 1];
        double x2 = mesh.V[n2 * 2], y2 = mesh.V[n2 * 2 + 1];
        /* Physical: (x,y) = (x0,y0) + xi*(x1-x0,y1-y0) + eta*(x2-x0,y2-y0) */
        for (int v = 0; v < n_sv; ++v) {
            double xi  = sv_xi[v], eta = sv_eta[v];
            double xp  = x0 + xi * (x1 - x0) + eta * (x2 - x0);
            double yp  = y0 + xi * (y1 - y0) + eta * (y2 - y0);
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
    std::cout << "[write_vtu] Solution written to: " << filepath << "\n";
}
