#include <stdexcept>
#include <string>
#include "quad.hpp"

typedef double real;

extern "C" {
	void shapeL(double *xref, int p, double **pphi);
	int gradientL(double *xref, int p, double **pgphi);

	extern int n1, n2, n4, n6;
	extern real x1[], x2[], x4[], x6[];
	extern real w1[], w2[], w4[], w6[];
}

QuadratureRule getQuadratureRule(int order) {
	switch (order) {
		case 0: return {n1, x1, w1};
		case 1: return {n2, x2, w2};
		case 2: return {n4, x4, w4};
		case 3: return {n6, x6, w6};
		default:
			throw std::runtime_error("Unsupported DG order: " + std::to_string(order));
	}
}

namespace {
	static const int n1d[] = {1, 2, 3, 4, 5};
	static const double x1d_1[] = {0.5};
	static const double w1d_1[] = {1.0};
	static const double x1d_2[] = {0.211324865405187, 0.788675134594813};
	static const double w1d_2[] = {0.5, 0.5};
	static const double x1d_3[] = {0.112701665379258, 0.5, 0.887298334620742};
	static const double w1d_3[] = {0.277777777777778, 0.444444444444444, 0.277777777777778};
	static const double x1d_4[] = {0.069431844202974, 0.330009478207572, 0.669990521792428, 0.930568155797026};
	static const double w1d_4[] = {0.173927422568727, 0.326072577431273, 0.326072577431273, 0.173927422568727};
	static const double x1d_5[] = {0.046910077030668, 0.230765344947158, 0.5, 0.769234655052841, 0.953089922969332};
	static const double w1d_5[] = {0.118463442528095, 0.239314335249683, 0.284444444444444, 0.239314335249683, 0.118463442528095};
}

QuadratureRule getQuadratureRule1D(int order) {
	if (order < 0 || order > 4)
		throw std::runtime_error("Unsupported 1D quadrature order: " + std::to_string(order));
	switch (order) {
		case 0: return {n1d[0], x1d_1, w1d_1};
		case 1: return {n1d[1], x1d_2, w1d_2};
		case 2: return {n1d[2], x1d_3, w1d_3};
		case 3: return {n1d[3], x1d_4, w1d_4};
		case 4: return {n1d[4], x1d_5, w1d_5};
		default: return {n1d[0], x1d_1, w1d_1};
	}
}

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
			int idx = i * quad.nq + q;
			phiq[idx] = {phi[i], gphi[i], gphi[i+Np], quad.wq[q]};
    }
	}
	if (phi) free(phi);
	if (gphi) free(gphi);
	return phiq;
}
