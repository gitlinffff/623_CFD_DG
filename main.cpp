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
    if (input.flux == "HLLE")
        flux_fn = fluxHLLE;
    else if (input.flux == "Roe")
        flux_fn = fluxROE;
    else if (input.flux == "Rusanov")
        flux_fn = fluxRusanov;
    else {
        std::cerr << "Unknown flux type\n";
        return 1;
    }

    const int Np = (order + 1) * (order + 2) / 2;
    std::vector<double> U(4 * mesh.Ne * Np, 0.0);
    initialize_uniform(U.data(), mesh.Ne, order, params);

    /* ---- Solve to steady state ---- */
    solve(mesh, U.data(), order, params, flux_fn,
                 CFL, print_interval, max_iter, use_local_dt);

    /* ---- Write converged solution to ParaView VTU ---- */
    std::string case_name = stem(gri_file) + "_p" + std::to_string(order);
    std::string outpath   = "results/steady/" + case_name + "/solution.vtu";
    write_solution_vtu(mesh, U.data(), order, params, outpath);
    std::cout << "Solution written to: " << outpath << "\n";

    return 0;
}
