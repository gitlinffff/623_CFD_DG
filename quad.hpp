#ifndef QUAD_HPP
#define QUAD_HPP

#include <vector>

// Structure to group quadrature data for easy passing
struct QuadratureRule {
	int nq;
	const double* xq;
	const double* wq;
};

struct BasisEval {
	double phi;
	double dphi_dxi;
	double dphi_deta;
	double wq;
};

/**
 * Returns the Dunavant quadrature rule appropriate for the given polynomial order.
 */
QuadratureRule getQuadratureRule(int order);

/**
 * Pre-computes basis functions and their gradients at all quadrature points.
 * Returns a vector with layout: [Basis i][Quad q][phi, dphi/dxi, dphi/deta, weight]
 */
std::vector<BasisEval> computePhiQ(int order);

#endif // QUAD_HPP
