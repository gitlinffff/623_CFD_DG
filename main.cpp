/**
 * Euler finite-volume solver - steady-state only.
 * Build: mkdir build && cd build && cmake .. && make
 * Usage: ./main
 * Output: data/results/solution.vtu
 */
#include "solver.hpp"
#include <iostream>
#include <vector>
#include <cstdio>
#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

static void ensure_dir(const char* path) {
#ifdef _WIN32
    (void)_mkdir(path);
#else
    (void)mkdir(path, 0755);
#endif
}

int run_steady(){
	int order = 2;  /* DG order*/
	int Np = (order + 1) * (order + 2) / 2;
	double Mach = 0.1;
	FluxFn flux_fn = fluxROE;
  ProblemParams params;

	const double t_end = 400;         /* run until periodic; adjust as needed */
	const double vtu_interval = 0.2;	
	const double CFL = 0.05;
	const char* gri_file = "/home/linfel/umich_course/623_CFD/mesh/ver2/global_refine_3.gri";
	const char* rst_file = "/home/linfel/umich_course/623_CFD/data/steady_results/refine3_2nd.vtu";
	const char* out_dir =  "/home/linfel/umich_course/623_CFD/data/unsteady/1st_unsteady_rfn3_solutions";

	GriMesh mesh;
	if (!read_gri(gri_file, mesh)) {
			std::cerr << "Failed to read mesh.\n";
			return 1;
	}

	std::vector<double> U(4 * mesh.Ne * Np);
	std::vector<double> R(4 * mesh.Ne * Np);
	initialize_uniform(U.data(), mesh.Ne, order, Mach, params);
  
	// this unfinished
	solve_steady(mesh, U.data(), gammad, params, flux_fn, recon_fn, CFL, 50, 1000000);
	return 0;
}


int rst_unsteady() {
    /* restart unsteady run*/
    const double gammad = 1.4;
    FluxFn flux_fn = fluxROE;
    ReconFn recon_fn = reconstruct_nolimiter;

    const double t_end = 400;         /* run until periodic; adjust as needed */
    const double vtu_interval = 0.2;	
    const double CFL = 0.05;
    const char* gri_file = "/home/linfel/umich_course/623_CFD/mesh/ver2/global_refine_3.gri";
    const char* rst_file = "/home/linfel/umich_course/623_CFD/data/steady_results/refine3_2nd.vtu";
    const char* out_dir =  "/home/linfel/umich_course/623_CFD/data/unsteady/1st_unsteady_rfn3_solutions";

    GriMesh mesh;
    if (!read_gri(gri_file, mesh)) {
        std::cerr << "Failed to read mesh.\n";
        return 1;
    }

    ensure_dir(out_dir);

    ProblemParams params;
    std::vector<double> U(mesh.Ne * 4);

    std::cout << "Mesh: " << mesh.Ne << " elements, "
              << mesh.num_interior_faces << " interior, "
              << mesh.num_boundary_faces << " boundary faces.\n";

    if (!read_vtu(mesh, rst_file, gammad, U.data())) {
        std::cerr << "Error: failed to read or mesh mismatch: " << rst_file << "\n";
        return 1;
    }
    std::cout << "Loaded: " << rst_file << "\n";

    solve_unsteady(mesh, U.data(), gammad, params, flux_fn, recon_fn, CFL, t_end, vtu_interval, 50, out_dir);

    return 0;
}

int run_unsteady() {
    /* run steady then unsteady based on steady solution */
    const double gammad = 1.4;
    FluxFn flux_fn = fluxROE;
    ReconFn recon_fn = reconstruct_nolimiter;

    const char* gri_file = "/home/linfel/umich_course/623_CFD/mesh/ver2/global_refine_1.gri";
    const char* out_dir =  "/home/linfel/umich_course/623_CFD/data/unsteady/test1";
    const double t_end = 300;         /* run until periodic; adjust as needed */
    const double vtu_interval = 0.2;	
    const double CFL = 0.3;

    GriMesh mesh;
    if (!read_gri(gri_file, mesh)) {
        std::cerr << "Failed to read mesh.\n";
        return 1;
    }
    ensure_dir(out_dir);

    ProblemParams params;
    std::vector<double> U(mesh.Ne * 4);

    std::cout << "Mesh: " << mesh.Ne << " elements, "
              << mesh.num_interior_faces << " interior, "
              << mesh.num_boundary_faces << " boundary faces.\n";

    initialize_uniform(U.data(), mesh.Ne, 0.1, params);

    solve_steady(mesh, U.data(), gammad, params, flux_fn, recon_fn, CFL, 50, 1000000);

    std::string filename = std::string(out_dir) + "/steady_solution.vtu";
    if (write_vtu(mesh, U.data(), gammad, filename.c_str()))
        std::cout << "Output: " << filename << "\n";
    else
        std::cerr << "Failed to write VTU.\n";

    solve_unsteady(mesh, U.data(), gammad, params, flux_fn, recon_fn, CFL, t_end, vtu_interval, 50, out_dir);
    
		return 0;
}
initialize_uniform_dg(U.data(), mesh.Ne, order, Mach, params);
int main() {
    /* switch run as needed */
    run_steady();
		//rst_unsteady();	
    //run_unsteady();
}
