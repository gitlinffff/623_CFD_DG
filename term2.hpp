#ifndef TERM2_HPP
#define TERM2_HPP

#include <vector>
#include <Eigen/Dense>
#include "physics.hpp"   // For ProblemParams
#include "quad.hpp"          // For QuadratureRule (assuming this defines the struct)
#include "readgri.hpp"
#include "problem.hpp"

/**
 * Computes the volume integral term (Term 2) for the DG residual:
 * Integral over element K of (Grad Phi_i \cdot PhysicalFlux)
 */
void addTerm2(const GridMesh& mesh, 
              std::vector<double>& R, 
              int order,
              const std::vector<double>& PhiQ, 
              const double* U,
              const ProblemParams& params);

#endif // TERM2_HPP
