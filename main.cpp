/**
 * DG solver main for the turbine blade passage.
 *
 * Edit the case parameters below and rebuild:
 *   order    : DG polynomial order (0, 1, 2, 3)
 *   CFL      : CFL number (recommend 0.4 for all orders)
 *   gri_file : path to .gri mesh file
 */
#include "advance.hpp"
#include "readgri.hpp"
#include "problem.hpp"
#include "write_vtu.hpp"
#include "parseinput.hpp"
#include "restart.hpp"
#include <algorithm>
#include <cctype>
#include <iostream>
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

    /* ---- Load mesh ---- */
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

    /* ---- Initialize solution to uniform freestream ---- */
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
        std::string restart_err;
        if (!initialize_from_restart(restart_data, order, mesh.Ne, U.data(), restart_err)) {
            std::cerr << "Failed to initialize from restart: " << restart_err << "\n";
            return 1;
        }
        std::cout << "Initialized from restart file: " << input.restart_file
                  << " (p" << restart_data.order << " -> p" << order << ")\n";
    } else {
        initialize_uniform(U.data(), mesh.Ne, order, params);
        std::cout << "Initialized from uniform freestream.\n";
    }

    std::string case_name = stem(gri_file) + "_p" + std::to_string(order);
    std::string case_dir  = "results/steady/" + case_name;
    std::string out_res   = case_dir + "/residual_history.dat";

    /* ---- Solve to steady state ---- */
    solve(mesh, U.data(), order, params, flux_fn,
                 CFL, print_interval, max_iter, use_local_dt, out_res);
    std::cout << "Residual history written to: " << out_res << "\n";

    /* ---- Write converged solution and restart coefficients ---- */
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
