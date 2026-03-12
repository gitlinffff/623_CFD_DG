#include "problem.hpp"
#include "physics.hpp"
#include <cstring>

namespace {
	const double PI = 3.14159265358979323846;
}

ProblemParams::ProblemParams() {
	rho0 = 1.0;
	a0 = 1.0;
	gammad = 1.4;
	alpha = 50.0 * PI / 180.0;  /* 50 deg */
	double p0 = getp0(rho0, a0, gammad);
	pout = 0.7 * p0;
	/* Unsteady wake (proj.pdf): Vrot=a0, delta_y=18mm, fwake=0.1, delta=0.1 */
	Vrot = 1.0;
	delta_y = 18;
	fwake = 0.1;
	delta_wake = 0.1;
}

double getp0(const ProblemParams& p) {
	return getp0(p.rho0, p.a0, p.gammad);
}


void initialize_uniform(double* U, int Ne, int order, const ProblemParams& params) {
	const double Mach = 0.1;  // freestream Mach number used for steady initialization

	int Np = (order + 1) * (order + 2) / 2;
	double rho, u, v, p;
	double p0 = getp0(params.rho0, params.a0, params.gammad);
	isentropic_prim_from_M(params.rho0, p0, params.gammad, Mach, params.alpha,
												 rho, u, v, p);

	// Convert to conserved variables
	double Ucons[4];
	primToCons(rho, u, v, p, params.gammad, Ucons);

	// Layout: [Element k][Variable var][Basis j]
	for (int var = 0; var < 4; ++var) {
		for (int k = 0; k < Ne; ++k) {
			for (int i = 0; i < Np; ++i) {
				int index = (var * Ne * Np) + (k * Np) + i;
				U[index] = Ucons[var];
			}
		}
	}
}
