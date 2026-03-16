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
						 double* dt_local, // Size: mesh.Ne
						 double& dt_global)
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

	// Boundary faces: same BC interpretation as FV solver
	addBndSurfTerm(mesh, R, order, U, params, flux_fn, sum_s);

	// Calculate dt_local and dt_global
	dt_global = 1.e20; // Initialize with large value
	#pragma omp parallel for schedule(static) reduction(min:dt_global)
	for (int k = 0; k < mesh.Ne; ++k) {
		double Ak = mesh.Area[k];

		double denominator = sum_s[k] * (2 * order + 1);
		
		if (denominator > 1e-15) {
			dt_local[k] = (CFL * 2.0 * Ak) / sum_s[k] / (2 * order + 1);
			// The (2p + 1) factor is for DG stability
		} else {
			dt_local[k] = 1e-6; // Safety fallback
		}
		if (dt_local[k] < dt_global) {dt_global = dt_local[k];}
	}
	if (dt_global >= 1.e20) {
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
					 bool use_local_dt)
{
    const int Np = (order + 1) * (order + 2) / 2;
    const double gammad = params.gammad;

    std::vector<double> R(4 * mesh.Ne * Np);
    std::vector<double> U1(4 * mesh.Ne * Np);
    std::vector<double> U2(4 * mesh.Ne * Np);
		std::vector<double> dt_local(mesh.Ne);
		std::vector<double> dt_dummy(mesh.Ne);
		double dt_global, dt_global_dummy;
		std::cout << "Use local time stepping: " << use_local_dt << std::endl;

    // Initial residual (M^{-1} applied so the norm is in DOF-update units)
    calcRes(mesh, U, R.data(), order, params, flux_fn, CFL, dt_dummy.data(), dt_global_dummy);
    applyInverseMassMatrix(mesh, R.data(), order);
    double R0 = residual_L1_norm(mesh, R.data(), order);
    std::cout << "Initial DG L1 residual: " << R0 << std::endl;

    double t = 0.0;
    int step = 0;

    while (step < max_iter) {
        // SSP-RK3 stage 1: U1_k = U_k - dt_k * M_k^{-1} R_k
        calcRes(mesh, U, R.data(), order, params, flux_fn, CFL, dt_local.data(), dt_global);
        applyInverseMassMatrix(mesh, R.data(), order);

				// an internal function to determine if dt_local or dt_global to use
				auto get_dt = [&](int k) {
            return use_local_dt ? dt_local[k] : dt_global;
        };

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
        calcRes(mesh, U1.data(), R.data(), order, params, flux_fn, CFL, dt_local.data(), dt_global);
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
        calcRes(mesh, U2.data(), R.data(), order, params, flux_fn, CFL, dt_local.data(), dt_global);
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

        t += dt_global; // t represents physical time for global time stepping
        ++step;

        if (residual_stride > 0 && step % residual_stride == 0) {
            calcRes(mesh, U, R.data(), order, params, flux_fn, CFL, dt_dummy.data(), dt_global_dummy);
            applyInverseMassMatrix(mesh, R.data(), order);
            double R1 = residual_L1_norm(mesh, R.data(), order);
            std::cout << "Step " << step << "  t=" << t << "  L1=" << R1;
            if (R0 > 1e-30)
                std::cout << "  ratio=" << (R1 / R0);
            std::cout << std::endl;
            if (R1 < R0 * 1e-5) {
                std::cout << "DG steady converged: L1 < 1e-5 * R0.\n";
                break;
            }
        }
    }
}
