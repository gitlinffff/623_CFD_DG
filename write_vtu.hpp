#ifndef WRITE_VTU_HPP
#define WRITE_VTU_HPP

#include "readgri.hpp"
#include "problem.hpp"
#include <string>

/**
 * Write the converged DG solution to a ParaView VTK Unstructured Grid (.vtu) file.
 *
 * Visualization strategy:
 *   For each element the reference triangle is uniformly subdivided into
 *   (2*order)^2 sub-triangles for p>=1.
 *   For p=0, subdivision is clamped to 1 to keep a valid output cell.
 *   The DG polynomial is evaluated at each sub-vertex so that intra-element
 *   variation is visible in ParaView for p >= 1.
 *   Each element owns its own sub-nodes (shared nodes between elements are NOT
 *   merged), which correctly represents the discontinuous nature of the solution.
 *
 * Output fields (PointData): rho, u, v, p, Mach, entropy (s/s_ref), wall_marker.
 *
 * The output directory is created automatically.
 *
 * @param mesh     Mesh read by read_gri().
 * @param U        DG DOF array, layout U[var * Ne * Np + k * Np + i].
 * @param order    Polynomial order (0,1,2,3).
 * @param params   Problem parameters (for entropy reference).
 * @param filepath Full path including filename, e.g.
 *                 "results/steady/global_refine_0_p0/solution.vtu".
 */
void write_solution_vtu(const GriMesh& mesh, const double* U, int order,
                        const ProblemParams& params, const std::string& filepath,
                        bool log_message = true);

void write_pvd_manifest(const std::string& case_dir, const std::vector<double>& times);

#endif
