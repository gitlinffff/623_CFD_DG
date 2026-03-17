#include "advance.hpp"
#include "readgri.hpp"
#include "problem.hpp"
#include "write_vtu.hpp"
#include "parseinput.hpp"
#include "restart.hpp"
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <iostream>
#include <limits.h>
#include <string>
#include <vector>

static std::string stem(const std::string& path)
{
    size_t slash = path.find_last_of("/\\");
    std::string base = (slash == std::string::npos) ? path : path.substr(slash + 1);
    size_t dot = base.rfind('.');
    return (dot == std::string::npos) ? base : base.substr(0, dot);
}

static std::string upper_copy(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    return s;
}

static bool canonical_path(const std::string& path, std::string& out)
{
    char buf[PATH_MAX];
    if (!realpath(path.c_str(), buf)) return false;
    out = buf;
    return true;
}

static bool same_mesh_file(const std::string& a, const std::string& b)
{
    if (a.empty() || b.empty()) return false;
    std::string ca, cb;
    if (canonical_path(a, ca) && canonical_path(b, cb)) return ca == cb;
    return a == b;
}

int main() {
    InputParams input;

    if (!read_input_file("input.dat", input))
        return 1;

    const int order = input.order;
    const double CFL = input.CFL;
    const char* gri_file = input.gri_file.c_str();
    const int max_iter = input.max_iter;
    const int print_interval = input.print_interval;
    const bool use_local_dt = input.use_local_dt;

    GriMesh mesh;
    if (!read_gri(gri_file, mesh)) {
        std::cerr << "Failed to read mesh: " << gri_file << "\n";
        return 1;
    }
    std::cout << "Mesh: " << gri_file
              << "  Ne=" << mesh.Ne
              << "  Nf_int=" << mesh.num_interior_faces
              << "  Nf_bnd=" << mesh.num_boundary_faces << "\n";
    std::cout << "Order p=" << order
              << "  CFL=" << CFL
							<< "  Flux=" << input.flux
							<< "  local_dt or not = " << use_local_dt << "\n";

    ProblemParams params;

    FluxFn flux_fn;
    const std::string flux_key = upper_copy(input.flux);
    if (flux_key == "HLLE")
        flux_fn = fluxHLLE;
    else if (flux_key == "HLLC")
        flux_fn = fluxHLLC;
    else if (flux_key == "ROE")
        flux_fn = fluxROE;
    else if (flux_key == "RUSANOV")
        flux_fn = fluxRusanov;
    else {
        std::cerr << "Unknown flux type: " << input.flux
                  << ". Supported: HLLE, HLLC, Roe, Rusanov\n";
        return 1;
    }

    const int Np = (order + 1) * (order + 2) / 2;
    std::vector<double> U(4 * mesh.Ne * Np, 0.0);
    if (!input.restart_file.empty()) {
        RestartData restart_data;
        if (!read_restart_dat(input.restart_file, restart_data)) {
            std::cerr << "Failed to read restart file: " << input.restart_file << "\n";
            return 1;
        }

        const bool same_mesh =
            same_mesh_file(restart_data.mesh_file, input.gri_file) ||
            (restart_data.mesh_file.empty() && restart_data.Ne == mesh.Ne);

        std::string restart_err;
        if (same_mesh) {
            if (!(order == restart_data.order || order == restart_data.order + 1)) {
                std::cerr << "Invalid same-mesh restart: only p->p or p->p+1 is allowed "
                          << "(got p" << restart_data.order << " -> p" << order << ").\n";
                return 1;
            }
            if (!initialize_from_restart(restart_data, order, mesh.Ne, U.data(), restart_err)) {
                std::cerr << "Failed to initialize from restart: " << restart_err << "\n";
                return 1;
            }
            std::cout << "Initialized from same-mesh restart: " << input.restart_file
                      << " (p" << restart_data.order << " -> p" << order << ")\n";
        } else {
            if (order != restart_data.order) {
                std::cerr << "Invalid cross-mesh restart: order must match exactly "
                          << "(got p" << restart_data.order << " -> p" << order << ").\n";
                return 1;
            }
            if (restart_data.mesh_file.empty()) {
                std::cerr << "Invalid cross-mesh restart: restart file does not contain mesh_file.\n";
                return 1;
            }
            GriMesh src_mesh;
            if (!read_gri(restart_data.mesh_file.c_str(), src_mesh)) {
                std::cerr << "Failed to read source mesh from restart metadata: "
                          << restart_data.mesh_file << "\n";
                return 1;
            }
            if (!initialize_from_restart_cross_mesh(restart_data, src_mesh, order, mesh,
                                                    U.data(), restart_err)) {
                std::cerr << "Failed to initialize from cross-mesh restart: " << restart_err << "\n";
                return 1;
            }
            std::cout << "Initialized from cross-mesh restart: " << input.restart_file
                      << " (mesh " << restart_data.mesh_file << " -> " << input.gri_file
                      << ", p" << restart_data.order << ")\n";
        }
    } else {
        initialize_uniform(U.data(), mesh.Ne, order, params);
        std::cout << "Initialized from uniform freestream.\n";
    }

    std::string case_name = stem(gri_file) + "_p" + std::to_string(order);
    std::string case_dir  = "results/steady/" + case_name;
    std::string out_res   = case_dir + "/residual_history.dat";

    solve(mesh, U.data(), order, params, flux_fn,
                 CFL, print_interval, max_iter, use_local_dt, out_res);
    std::cout << "Residual history written to: " << out_res << "\n";

    std::string out_vtu   = case_dir + "/solution.vtu";
    std::string out_dat   = case_dir + "/restart.dat";

    write_solution_vtu(mesh, U.data(), order, params, out_vtu);
    std::cout << "Solution written to: " << out_vtu << "\n";
    if (!write_restart_dat(out_dat, mesh, order, U.data(), input.gri_file)) {
        std::cerr << "Failed to write restart file: " << out_dat << "\n";
        return 1;
    }
    std::cout << "Restart coefficients written to: " << out_dat << "\n";

    return 0;
}
