#ifndef ADVANCE_HPP
#define ADVANCE_HPP

#include "readgri.hpp"
#include "problem.hpp"
#include "reconstruct.hpp"

/** Flux function: flux(UL, UR, n, gammad, Fhat, smag) */
typedef void (*FluxFn)(const double*, const double*, const double*, double, double*, double&);

/**
 * Compute semi-discrete residual R = -(1/A_i) * sum_faces (Fhat * L).
 * Interior: flux(UL, UR, n). Boundary: WallFlux / InflowFlux / OutflowFlux per type.
 * recon_fn: reconstruction (e.g. reconstruct_const).
 * If dt_per_cell != nullptr, also compute local dt_i = (2*A_i*CFL)/(sum_e |s|_e*L_e) per cell.
 */
void calcRes(const GriMesh& mesh, const double* U, double* R, double gammad,
            const ProblemParams& params, FluxFn flux_fn, ReconFn recon_fn,
            double* dt_per_cell, double& dt_min, double CFL = 0.5, double t = -1.0);

/** SSP-RK3 with global dt (generic). Returns dt used. t: current time for time-dependent inflow. */
double SSPRK3(const GriMesh& mesh, double* U, double gammad,
              const ProblemParams& params, FluxFn flux_fn, ReconFn recon_fn,
              double CFL, double t);

/** SSP-RK3 with local time stepping (for steady-state convergence). */
void SSPRK3_local(const GriMesh& mesh, double* U, double gammad,
                  const ProblemParams& params, FluxFn flux_fn, ReconFn recon_fn,
                  double CFL = 0.5);

/** Unsteady solve: global dt, time-dependent inflow. Outputs to out_dir. */
void solve_unsteady(const GriMesh& mesh, double* U, double gammad,
                    const ProblemParams& params, FluxFn flux_fn, ReconFn recon_fn,
                    double CFL, double t_end, double vtu_interval, int residual_stride,
                    const char* out_dir);

double residual_L1_norm(const GriMesh& mesh, const double* R);

/** L2 norm of residual: sqrt(sum_i R_i^2) */
double residual_L2_norm(const GriMesh& mesh, const double* R);

/** Global CFL-based dt (for reference). Local time stepping uses dt_per_cell from calcRes. */
double compute_dt(const GriMesh& mesh, const double* U, double gammad, double CFL);

/**
 * Steady-state solve: SSP-RK3 until L1 residual drops 5 orders below initial.
 * flux_fn: fluxROE, fluxHLLC, or fluxRusanov
 * recon_fn: reconstruct_const or reconstruct_nolimiter
 */
void solve_steady(const GriMesh& mesh, double* U, double gammad, const ProblemParams& params,
                  FluxFn flux_fn, ReconFn recon_fn, double CFL = 0.1, int residual_stride = 50,
                  int max_iter = 100000);

#endif
