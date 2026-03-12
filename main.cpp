/**
 * DG solver main for the turbine blade passage.
 *
 * Usage: build with 'make' from the project root, run with 'make run'.
 *
 * Edit the three parameters below to change the case:
 *   order    : DG polynomial order (0 = FV-equivalent, 1, 2, 3)
 *   gri_file : path to .gri mesh file (relative to project root)
 */
#include "advance.hpp"
#include "readgri.hpp"
#include "problem.hpp"
#include "write_vtu.hpp"
#include <iostream>
#include <string>
#include <vector>

/* Extract the stem of a file path (no directory prefix, no extension). */
static std::string stem(const std::string& path)
{
    size_t slash = path.find_last_of("/\\");
    std::string base = (slash == std::string::npos) ? path : path.substr(slash + 1);
    size_t dot = base.rfind('.');
    return (dot == std::string::npos) ? base : base.substr(0, dot);
}

int main()
{
    /* ---- Case parameters ---- */
    const int    order    = 0;
    const double CFL      = 0.4;
    const char*  gri_file = "mesh/global_refine_0.gri";

    /* ---- Load mesh ---- */
    GriMesh mesh;
    if (!read_gri(gri_file, mesh)) {
        std::cerr << "Failed to read mesh: " << gri_file << "\n";
        return 1;
    }
    std::cout << "Mesh: " << gri_file
              << "  (" << mesh.Ne << " elements, "
              << mesh.num_interior_faces << " interior faces)\n";

    /* ---- Initialise solution ---- */
    ProblemParams params;
    const int Np = (order + 1) * (order + 2) / 2;
    std::vector<double> U(4 * mesh.Ne * Np, 0.0);
    initialize_uniform(U.data(), mesh.Ne, order, params);

    /* ---- Steady solve ---- */
    FluxFn flux_fn = fluxROE;
    solve_steady(mesh, U.data(), order, params, flux_fn, CFL, /*print_interval=*/50, /*max_iter=*/100000);

    /* ---- Write converged solution to ParaView VTU ---- */
    std::string case_name = stem(gri_file) + "_p" + std::to_string(order);
    std::string outpath   = "results/steady/" + case_name + "/solution.vtu";
    write_solution_vtu(mesh, U.data(), order, params, outpath);

    return 0;
}
