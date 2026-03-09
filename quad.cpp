#include <stdexcept>
#include <string>
#include "quad.hpp"

// Define 'real' (as in quad2d.c)
typedef double real;

extern "C" {
	void shapeL(double *xref, int p, double **pphi);
	int gradientL(double *xref, int p, double **pgphi);

	// Dunavant variables from quad2d.c
	extern int n1, n2, n4, n6;
	extern real x1[], x2[], x4[], x6[];
	extern real w1[], w2[], w4[], w6[];
}

QuadratureRule getQuadratureRule(int order) {
	// The integration of phi_i * phi_j requires a rule of degree 2*p
	switch (order) {
		case 0: return {n1, x1, w1};
		case 1: return {n2, x2, w2};
		case 2: return {n4, x4, w4};
		case 3: return {n6, x6, w6};
		default:
			throw std::runtime_error("Unsupported DG order: " + std::to_string(order));
	}
}

/* get phi, dPhi/dXi, dPhi/dEta, wq at all quad points for all basis functions*/
std::vector<BasisEval> computePhiQ(int order) {
	int Np = (order + 1) * (order + 2) / 2;
	QuadratureRule quad = getQuadratureRule(order);
	std::vector<BasisEval> phiq(Np * quad.nq);
	double xref[2];
	double *phi = nullptr;
	double *gphi = nullptr;
	
	for (int q = 0; q < quad.nq; ++q) {
		xref[0] = quad.xq[2 * q];
		xref[1] = quad.xq[2 * q + 1];
		shapeL(xref, order, &phi);
		gradientL(xref, order, &gphi);
		for (int i = 0; i < Np; ++i) {
			// Layout: [Basis i][Quad q][BasisEval]
			int idx = i * quad.nq + q;
			// Phi, dPhi/dXi, dPhi/dEta, quad weight
			phiq[idx] = {phi[i], gphi[i], gphi[i+Np], quad.wq[q]};
    }
	}
	if (phi) free(phi);
	if (gphi) free(gphi);
	return phiq;
}
