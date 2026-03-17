/**
 * DG time advance and residual for the turbine-passage Euler solver.
 *
 * This replaces the original FV-based advance.cpp with a DG formulation:
 *  - Residuals are assembled with addTerm2/addSurfTerm/addBndSurfTerm.
 *  - Time stepping uses SSP-RK3 with a global CFL-based dt.
 */
#include "advance.hpp"
#include "mm.hpp"
#include <cmath>
#include <cstring>
#include <iostream>
#include <vector>

void calcRes(const GriMesh& mesh,
						 const double* U,
						 double* R,
						 int order,
						 const ProblemParams& params,
						 FluxFn flux_fn,
						 double CFL,
						 double* dt_loc, // Size: mesh.Ne
						 double& dt_glb,
						 bool in_ptb,
						 const double t)
{
	const int Np = (order + 1) * (order + 2) / 2;

	// Zero-out residual
	std::memset(R, 0, sizeof(double) * 4 * mesh.Ne * Np);
	std::vector<double> sum_s(mesh.Ne, 0.0);// sum(max_smag_e * L_e) for 3 edges of each elem

	// Pre-compute basis and gradients at quadrature points (cached per order)
	static std::vector<BasisEval> phiq_cache;
	static int cached_order = -1;
	if (cached_order != order) {
			phiq_cache = computePhiQ(order);
			cached_order = order;
	}

	// Volume contribution
	addTerm2(mesh, R, order, phiq_cache, U, params);

	// Interior faces
	addSurfTerm(mesh, R, order, U, params, flux_fn, sum_s);

	// Boundary faces
	addBndSurfTerm(mesh, R, order, U, params, flux_fn, sum_s, in_ptb, t);

	// Calculate dt_loc and dt_glb
	dt_glb = 1.e20; // Initialize with large value
	#pragma omp parallel for schedule(static) reduction(min:dt_glb)
	for (int k = 0; k < mesh.Ne; ++k) {
		double Ak = mesh.Area[k];

		double denominator = sum_s[k] * (2 * order + 1);
		
		if (denominator > 1e-15) {
			dt_loc[k] = (CFL * 2.0 * Ak) / sum_s[k] / (2 * order + 1);
			// The (2p + 1) factor is for DG stability
		} else {
			dt_loc[k] = 1e-6; // Safety fallback
		}
		if (dt_loc[k] < dt_glb) {dt_glb = dt_loc[k];}
	}
	if (dt_glb >= 1.e20) {
		throw std::runtime_error("Global minimum dt calculation failed.");
	}
}

double residual_L1_norm(const GriMesh& mesh,
                           const double* R,
                           int order)
{
    const int Np = (order + 1) * (order + 2) / 2;
    const int nTot = 4 * mesh.Ne * Np;
    double sum = 0.0;
    #pragma omp parallel for schedule(static) reduction(+:sum)
    for (int i = 0; i < nTot; ++i)
        sum += std::fabs(R[i]);
    return sum;
}

void solve(const GriMesh& mesh,
					 double* U,
					 int order,
					 const ProblemParams& params,
					 FluxFn flux_fn,
					 double CFL,
					 int residual_stride,
					 int max_iter,
					 bool use_local_dt,
					 bool in_ptb)
{
	if (in_ptb && use_local_dt) {
			throw std::runtime_error(
					"Must use global time stepping when inflow unsteady perturbation is on.");
	}

	const int Np = (order + 1) * (order + 2) / 2;
	const double gammad = params.gammad;

	std::vector<double> R(4 * mesh.Ne * Np);
	std::vector<double> U1(4 * mesh.Ne * Np);
	std::vector<double> U2(4 * mesh.Ne * Np);
	std::vector<double> dt_loc(mesh.Ne);
	std::vector<double> dt_dmy(mesh.Ne); // dt_dummy
	double dt_glb, dt_gdmy; // dt_global_dummy

	double R0;
	double t = 0.0;
	int step = 0;

	while (step < max_iter) {
		// calculate residuals at current state
		calcRes(mesh, U, R.data(), order, params, flux_fn, CFL, dt_loc.data(), dt_glb, in_ptb, t);
		
		if (residual_stride > 0 && step % residual_stride == 0) {
			double R1 = residual_L1_norm(mesh, R.data(), order);
			if (step == 0) R0 = R1;
			std::cout << "Step " << step << "  t=" << t << "  L1=" << R1;
			if (R0 > 1e-30)
				std::cout << "  ratio=" << (R1 / R0);
			std::cout << std::endl;
			if (R1 < R0 * 1e-5) {
				std::cout << "DG steady converged: L1 < 1e-5 * R0.\n";
				break;
			}
		}
	
		auto get_dt = [&](int k) {
			return use_local_dt ? dt_loc[k] : dt_glb;
		}; // an internal function to determine if dt_loc or dt_glb to use

		// SSP-RK3 stage 1: U1_k = U_k - dt_k * M_k^{-1} R_k
		applyInverseMassMatrix(mesh, R.data(), order);
		#pragma omp parallel for schedule(static)
		for (int k = 0; k < mesh.Ne; ++k) {
			double dtk = get_dt(k);
			for (int var = 0; var < 4; ++var) {
				int base = (var * mesh.Ne + k) * Np;
				for (int i = 0; i < Np; ++i)
					U1[base + i] = U[base + i] - dtk * R[base + i];
			}
		}

		// stage 2: U2_k = 0.75 U_k + 0.25 (U1_k - dt_k * M_k^{-1} R1_k)
		calcRes(mesh, U1.data(), R.data(), order, params, flux_fn, CFL, dt_dmy.data(), dt_gdmy, in_ptb, t+dt_glb);
		applyInverseMassMatrix(mesh, R.data(), order);
		#pragma omp parallel for schedule(static)
		for (int k = 0; k < mesh.Ne; ++k) {
			double dtk = get_dt(k);
			for (int var = 0; var < 4; ++var) {
				int base = (var * mesh.Ne + k) * Np;
				for (int i = 0; i < Np; ++i)
					U2[base + i] = 0.75 * U[base + i] + 0.25 * (U1[base + i] - dtk * R[base + i]);
			}
		}

		// stage 3: U_k = 1/3 U_k + 2/3 (U2_k - dt_k * M_k^{-1} R2_k)
		calcRes(mesh, U2.data(), R.data(), order, params, flux_fn, CFL, dt_dmy.data(), dt_gdmy, in_ptb, t+0.5*dt_glb);
		applyInverseMassMatrix(mesh, R.data(), order);
		#pragma omp parallel for schedule(static)
		for (int k = 0; k < mesh.Ne; ++k) {
			double dtk = get_dt(k);
			for (int var = 0; var < 4; ++var) {
				int base = (var * mesh.Ne + k) * Np;
				for (int i = 0; i < Np; ++i)
					U[base + i] = (1.0/3.0)*U[base + i] + (2.0/3.0)*(U2[base + i] - dtk * R[base + i]);
			}
		}

		t += dt_glb; // t represents physical time for global time stepping
		++step;
	}
}
