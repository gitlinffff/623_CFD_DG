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
 * Returns 1D Gauss-Legendre quadrature rule for edge integration.
 * order 0 -> 1 pt, 1 -> 2 pts, 2 -> 3 pts, 3 -> 4 pts (for degree ~2p on edge).
 */
QuadratureRule getQuadratureRule1D(int order);

/**
 * Pre-computes basis functions and their gradients at all quadrature points.
 * Returns a vector with layout: [Basis i][Quad q][phi, dphi/dxi, dphi/deta, weight]
 */
std::vector<BasisEval> computePhiQ(int order);

#endif // QUAD_HPP
