#ifndef SURF_HPP
#define SURF_HPP

#include "readgri.hpp"
#include "problem.hpp"
#include "quad.hpp"
#include "flux.hpp"

typedef void (*FluxFn)(const double*, const double*, const double*, double, double*, double&);

void addSurfTerm(const GriMesh& mesh, double* R, int order,
                const double* U, const ProblemParams& params, FluxFn flux_fn,
								std::vector<double>& sum_s);

enum class DGBcType {
    WALL,
    INFLOW,
    OUTFLOW,
    FREESTREAM
};

DGBcType classify_boundary_name(const std::string& name_raw);

void addBndSurfTerm(const GriMesh& mesh, double* R, int order,
                    const double* U, const ProblemParams& params, FluxFn flux_fn,
										std::vector<double>& sum_s, bool in_ptb, const double t);

#endif
