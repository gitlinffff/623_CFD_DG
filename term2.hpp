#ifndef TERM2_HPP
#define TERM2_HPP

#include <vector>
#include <Eigen/Dense>
#include "quad.hpp"
#include "readgri.hpp"
#include "problem.hpp"

/**
 * Computes the volume integral term (Term 2) for the DG residual:
 * Integral over element K of (Grad Phi_i \cdot PhysicalFlux)
 */
void addTerm2(const GriMesh& mesh, double* R, int order,
	            const std::vector<BasisEval>& phiq, const double* U,
							const ProblemParams& params);

#endif // TERM2_HPP
