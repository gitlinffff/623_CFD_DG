#ifndef TERM2_HPP
#define TERM2_HPP

#include <vector>
#include <Eigen/Dense>
#include "quad.hpp"
#include "readgri.hpp"
#include "problem.hpp"

void addTerm2(const GriMesh& mesh, double* R, int order,
	            const std::vector<BasisEval>& phiq, const double* U,
							const ProblemParams& params);

#endif
