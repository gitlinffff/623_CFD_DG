#ifndef ADVANCE_HPP
#define ADVANCE_HPP

#include "readgri.hpp"
#include "problem.hpp"
#include "flux.hpp"
#include "quad.hpp"
#include "term2.hpp"
#include "surf.hpp"
#include "physics.hpp"
#include <string>

typedef void (*FluxFn)(const double*, const double*, const double*, double, double*, double&);

void calcRes(const GriMesh& mesh,
             const double* U,
             double* R,
             int order,
             const ProblemParams& params,
             FluxFn flux_fn,
             double CFL,
             double* dt_local,
             double& dt_global);

double residual_L1_norm(const GriMesh& mesh,
                           const double* R,
                           int order);

void solve(const GriMesh& mesh,
                  double* U,
                  int order,
                  const ProblemParams& params,
                  FluxFn flux_fn,
                  double CFL = 0.5,
                  int residual_stride = 50,
                  int max_iter = 1000000,
									bool use_local_dt = false,
                  const std::string& residual_history_file = "");

#endif
