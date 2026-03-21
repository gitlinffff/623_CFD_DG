#ifndef SURF_HPP
#define SURF_HPP

#include "readgri.hpp"
#include "problem.hpp"
#include "quad.hpp"
#include "flux.hpp"

/** Flux function: flux(UL, UR, n, gammad, Fhat, smag) */
typedef void (*FluxFn)(const double*, const double*, const double*, double, double*, double&);

/**
 * Computes the surface flux integral term for the DG residual (interior faces only).
 * R_i += int_e phi_i * Fhat(U-, U+, n) ds  for elemL
 * R_i -= int_e phi_i * Fhat(U-, U+, n) ds  for elemR  (conservation: Fhat(UR,UL,-n) = -Fhat(UL,UR,n))
 *
 * Boundary faces are not handled (leave empty as requested).
 */
void addSurfTerm(const GriMesh& mesh, double* R, int order,
                const double* U, const ProblemParams& params, FluxFn flux_fn,
								std::vector<double>* sum_s = nullptr);

/**
 * Computes the surface flux integral term over boundary faces for the DG residual.
 *
 * Default behavior (freestream test): use a ghost state equal to the interior state,
 * i.e. Fhat = flux(UL, UL, n). This makes the constant freestream solution a
 * discrete steady state when combined with the volume term.
 */
void addBndSurfTerm(const GriMesh& mesh, double* R, int order,
                    const double* U, const ProblemParams& params, FluxFn flux_fn,
										bool in_ptb, const double t, std::vector<double>* sum_s = nullptr);

#endif // SURF_HPP
