#ifndef ADVANCE_HPP
#define ADVANCE_HPP

#include "readgri.hpp"
#include "problem.hpp"
#include "flux.hpp"
#include "quad.hpp"
#include "term2.hpp"
#include "surf.hpp"
#include "physics.hpp"

/** Flux function: flux(UL, UR, n, gammad, Fhat, smag) */
typedef void (*FluxFn)(const double*, const double*, const double*, double, double*, double&);

/**
 * Compute DG residual R(U) on the given mesh.
 *
 * - U layout: [var][elem][basis] with size 4 * Ne * Np,
 *   where Np = (order+1)(order+2)/2.
 * - R has the same layout and size.
 *
 * Uses addTerm2 (volume), addSurfTerm (interior), and addBndSurfTerm (boundary)
 * so that the DG solver sees the same geometry and BCs as the old FV solver.
 */
void calcRes(const GriMesh& mesh,
             const double* U,
             double* R,
             int order,
             const ProblemParams& params,
             FluxFn flux_fn,
             double CFL,
             double* dt_local, // Size: mesh.Ne
             double& dt_global);


/** L1 norm of DG residual over all variables/elements/basis functions. */
double residual_L1_norm(const GriMesh& mesh,
                           const double* R,
                           int order);

/**
 * Steady-state DG solve: SSP-RK3 with (currently) local time stepping.
 *
 * Marches in pseudo-time until the L1 norm of the residual drops 5 orders
 * of magnitude relative to the initial value, mirroring the FV solver logic.
 */
void solve(const GriMesh& mesh,
                  double* U,
                  int order,
                  const ProblemParams& params,
                  FluxFn flux_fn,
                  double CFL = 0.5,
                  int residual_stride = 50,
                  int max_iter = 1000000,
									bool use_local_dt = false);

#endif
