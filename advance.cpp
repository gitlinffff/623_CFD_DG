#include "advance.hpp"
#include "restart.hpp"
#include "write_vtu.hpp"
#include "mm.hpp"
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sys/stat.h>
#include <vector>

namespace {
void make_dirs(const std::string& path)
{
    for (size_t i = 1; i <= path.size(); ++i) {
        if (i == path.size() || path[i] == '/') {
            std::string sub = path.substr(0, i);
            mkdir(sub.c_str(), 0755);
        }
    }
}
}

void calcRes(const GriMesh& mesh,
						 const double* U,
						 double* R,
						 int order,
						 const ProblemParams& params,
						 FluxFn flux_fn,
						 double CFL,
						 double* dt_loc,
						 double& dt_glb,
						 bool in_ptb,
						 const double t,
						 const std::map<int, ElementMetrics>& curved_metrics)
{
	const int Np = (order + 1) * (order + 2) / 2;

	std::memset(R, 0, sizeof(double) * 4 * mesh.Ne * Np);
	std::vector<double> sum_s(mesh.Ne, 0.0);

	static std::vector<BasisEval> phiq_cache;
	static int cached_order = -1;
	if (cached_order != order) {
			phiq_cache = computePhiQ(order);
			cached_order = order;
	}

	addTerm2(mesh, R, order, phiq_cache, U, params, curved_metrics);

	addSurfTerm(mesh, R, order, U, params, flux_fn, sum_s);

	addBndSurfTerm(mesh, R, order, U, params, flux_fn, sum_s, in_ptb, t);

	dt_glb = 1.e20;
	#pragma omp parallel for schedule(static) reduction(min:dt_glb)
	for (int k = 0; k < mesh.Ne; ++k) {
		double Ak = mesh.Area[k];

		double denominator = sum_s[k] * (2 * order + 1);

		if (denominator > 1e-15) {
			dt_loc[k] = (CFL * 2.0 * Ak) / sum_s[k] / (2 * order + 1);
		} else {
			dt_loc[k] = 1e-6;
		}
		if (dt_loc[k] < dt_glb) {dt_glb = dt_loc[k];}
	}
	if (dt_glb >= 1.e20) {
		throw std::runtime_error("Global minimum dt calculation failed.");
	}
}

double residual_L1_norm(const GriMesh& mesh, const double* R, int order) {
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
					 bool in_ptb,
           const std::string& residual_history_file,
					 const std::string& case_dir,
					 const double vtu_interval,
					 const double dat_interval,
					 const double t_final,
					 const std::string& mesh_file)
{
	if (in_ptb && use_local_dt) {
		throw std::runtime_error(
			"Must use global time stepping when inflow unsteady perturbation is on.");
	}

	std::map<int, ElementMetrics> curved_metrics;
	std::vector<BasisEval> phiq_cache = computePhiQ(order);
	std::cout << "=========================" << std::endl;
	precomputeCurvedMetrics(mesh, order, phiq_cache, curved_metrics);

	const int Np = (order + 1) * (order + 2) / 2;
	const double gammad = params.gammad;

	std::vector<double> R(4 * mesh.Ne * Np);
	std::vector<double> U1(4 * mesh.Ne * Np);
	std::vector<double> U2(4 * mesh.Ne * Np);
	std::vector<double> dt_loc(mesh.Ne);
	std::vector<double> dt_dmy(mesh.Ne); 
	double dt_glb, dt_gdmy;

	double R0 = -1.0;
	double t = 0.0;
	int step = 0;
	double next_vtu_t = vtu_interval;
	double next_dat_t = dat_interval;
	std::vector<double> vtu_times;

	std::ofstream history_out;
	if (!residual_history_file.empty()) {
		size_t slash = residual_history_file.find_last_of('/');
		if (slash != std::string::npos)
			make_dirs(residual_history_file.substr(0, slash));

		history_out.open(residual_history_file.c_str());
		if (!history_out) {
			throw std::runtime_error("Failed to open residual history file: " + residual_history_file);
		}
		history_out << "# step t L1 ratio\n";
		history_out << std::setprecision(17);
	}

	while (step < max_iter && (!in_ptb || t < t_final - 1e-12)) {
		calcRes(mesh, U, R.data(), order, params, flux_fn, CFL,
		        dt_loc.data(), dt_glb, in_ptb, t, curved_metrics);

		if (residual_stride > 0 && step % residual_stride == 0) {
			double R1 = residual_L1_norm(mesh, R.data(), order);
			if (step == 0) R0 = R1;
			double ratio = (R0 > 1e-30) ? (R1 / R0) : 1.0;

			if (history_out)
				history_out << step << " " << t << " " << R1 << " " << ratio << "\n";

			std::cout << "Step " << step << "  t=" << t << "  L1=" << R1;
			if (R0 > 1e-30)
				std::cout << "  ratio=" << ratio;
			std::cout << std::endl;
			if (R1 < R0 * 1e-5) {
				std::cout << "DG steady converged: L1 < 1e-5 * R0.\n";
				break;
			}
		}

		auto get_dt = [&](int k) {
			return use_local_dt ? dt_loc[k] : dt_glb;
		};

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

		calcRes(mesh, U1.data(), R.data(), order, params, flux_fn, CFL, dt_dmy.data(),
		        dt_gdmy, in_ptb, t+dt_glb, curved_metrics);
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

		calcRes(mesh, U2.data(), R.data(), order, params, flux_fn, CFL, dt_dmy.data(),
		        dt_gdmy, in_ptb, t+0.5*dt_glb, curved_metrics);
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

		t += dt_glb;
		++step;

		if (in_ptb) {
			if (dat_interval > 0. && t >= next_dat_t - 1e-12) { // .dat output
					char dat_path[256];
					std::snprintf(dat_path, sizeof(dat_path), "%s/restart_t_%05.0f.dat", case_dir.c_str(), t);
					write_restart_dat(dat_path, mesh, order, U, mesh_file);
					next_dat_t += dat_interval;
			}

			if (vtu_interval > 0. && t >= next_vtu_t - 1e-12) { // .vtu output
					char vtu_path[256];
					std::snprintf(vtu_path, sizeof(vtu_path), "%s/solution_t_%05.0f.vtu", case_dir.c_str(), t);
					write_solution_vtu(mesh, U, order, params, vtu_path);
					vtu_times.push_back(t);
					write_pvd_manifest(case_dir, vtu_times);
					next_vtu_t += vtu_interval;
			}
		}
	}
	
}
