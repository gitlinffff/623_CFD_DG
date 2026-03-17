#ifndef QUAD_HPP
#define QUAD_HPP

#include <vector>

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

QuadratureRule getQuadratureRule(int order);

QuadratureRule getQuadratureRule1D(int order);

std::vector<BasisEval> computePhiQ(int order);

#endif
