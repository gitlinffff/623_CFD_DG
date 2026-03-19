#ifndef TERM2_HPP
#define TERM2_HPP

#include <vector>
#include <Eigen/Dense>
#include <map>
#include "quad.hpp"
#include "readgri.hpp"
#include "problem.hpp"


struct ElementMetrics {
    std::vector<Eigen::Matrix2d> G; // store J_inv * detJ
    std::vector<double> detJ;
};

void precomputeCurvedMetrics(const GriMesh& mesh, int order,
                            const std::vector<BasisEval>& phiq,
                            std::map<int, ElementMetrics>& curved_metrics);

void addTerm2(const GriMesh& mesh, double* R, int order,
              const std::vector<BasisEval>& phiq, const double* U,
              const ProblemParams& params,
              const std::map<int, ElementMetrics>& curved_metrics);

#endif
