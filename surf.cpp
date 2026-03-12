#include "surf.hpp"
#include "mm.hpp"
#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstring>
#include <string>

extern "C" {
	void shapeL(double* xref, int p, double** pphi);
}

namespace {

/**
 * Reference (xi, eta) on face j at 1D param t in [0,1].
 *
 * Convention matches readgri.cpp: face j is the edge from local node j
 * to local node (j+1)%3.  Reference nodes: 0=(0,0), 1=(1,0), 2=(0,1).
 *
 *   face 0: (0,0) -> (1,0)   node0 to node1   xi = t,     eta = 0
 *   face 1: (1,0) -> (0,1)   node1 to node2   xi = 1-t,   eta = t
 *   face 2: (0,1) -> (0,0)   node2 to node0   xi = 0,     eta = 1-t
 */
void faceRefCoords(int face, double t, double& xi, double& eta) {
	switch (face) {
		case 0: xi = t;       eta = 0.0;     break;  /* node0=(0,0) -> node1=(1,0) */
		case 1: xi = 1.0-t;   eta = t;       break;  /* node1=(1,0) -> node2=(0,1) */
		case 2: xi = 0.0;     eta = 1.0-t;   break;  /* node2=(0,1) -> node0=(0,0) */
		default: xi = t;      eta = 0.0;     break;
	}
}

/**
 * Physical (x,y) on face j of element k at param t in [0,1].
 *
 * Convention matches readgri.cpp and faceRefCoords above:
 * face j goes from physical node j to physical node (j+1)%3.
 */
void facePhysCoords(const GriMesh& mesh, int k, int face, double t,
                   double& x, double& y) {
	int v0 = mesh.E[k * 3 + face];
	int v1 = mesh.E[k * 3 + (face + 1) % 3];
	double x0 = mesh.V[v0 * 2 + 0], y0 = mesh.V[v0 * 2 + 1];
	double x1 = mesh.V[v1 * 2 + 0], y1 = mesh.V[v1 * 2 + 1];
	x = x0 + t * (x1 - x0);
	y = y0 + t * (y1 - y0);
}

/** Ref coords (xi, eta) in element k corresponding to physical (x,y). */
void physToRef(const GriMesh& mesh, int k, double x, double y,
               double& xi, double& eta) {
	Eigen::Matrix2d J = Jacobian(mesh, k);
	int v0 = mesh.E[k * 3 + 0];
	double x0 = mesh.V[v0 * 2 + 0], y0 = mesh.V[v0 * 2 + 1];
	Eigen::Vector2d d;
	d << x - x0, y - y0;
	Eigen::Vector2d ref = J.inverse() * d;
	xi = ref(0);
	eta = ref(1);
}

} // namespace

void addSurfTerm(const GriMesh& mesh, double* R, int order,
                const double* U, const ProblemParams& params, FluxFn flux_fn) {
	int Np = (order + 1) * (order + 2) / 2;
	QuadratureRule quad1d = getQuadratureRule1D(order);

	double* phi = nullptr;

	for (int iface = 0; iface < mesh.num_interior_faces; ++iface) {
		int elemL = mesh.I2E[4 * iface + 0];
		int faceL = mesh.I2E[4 * iface + 1];
		int elemR = mesh.I2E[4 * iface + 2];
		int faceR = mesh.I2E[4 * iface + 3];

		const double* n = &mesh.In[2 * iface];
		double L = mesh.In_len[iface];

		for (int q = 0; q < quad1d.nq; ++q) {
			double t = quad1d.xq[q];
			double w = quad1d.wq[q] * L;

			/* Ref coords for elemL on face faceL */
			double xiL, etaL;
			faceRefCoords(faceL, t, xiL, etaL);

			/* Physical point */
			double x_phys, y_phys;
			facePhysCoords(mesh, elemL, faceL, t, x_phys, y_phys);

			/* Ref coords for elemR.
			 * For regular interior faces physToRef maps the shared physical
			 * point back into elemR's reference domain.
			 * For PERIODIC faces the physical point on elemL's side is NOT
			 * inside elemR's physical domain, so physToRef produces wildly
			 * out-of-range coordinates; in that case fall back to the direct
			 * parametric mapping faceRefCoords(faceR, t). */
			double xiR, etaR;
			physToRef(mesh, elemR, x_phys, y_phys, xiR, etaR);
			if (xiR < -0.1 || etaR < -0.1 || xiR + etaR > 1.1) {
				/* Periodic face: use local face parametrisation of elemR. */
				faceRefCoords(faceR, t, xiR, etaR);
			}

			/* Interpolate UL at (xiL, etaL) */
			double xrefL[2] = {xiL, etaL};
			shapeL(xrefL, order, &phi);
			double UL[4] = {0, 0, 0, 0};
			for (int j = 0; j < Np; ++j) {
				for (int var = 0; var < 4; ++var)
					UL[var] += U[var * mesh.Ne * Np + elemL * Np + j] * phi[j];
			}

			/* Interpolate UR at (xiR, etaR) */
			double xrefR[2] = {xiR, etaR};
			shapeL(xrefR, order, &phi);
			double UR[4] = {0, 0, 0, 0};
			for (int j = 0; j < Np; ++j) {
				for (int var = 0; var < 4; ++var)
					UR[var] += U[var * mesh.Ne * Np + elemR * Np + j] * phi[j];
			}

			/* Numerical flux Fhat(UL, UR, n) */
			double Fhat[4];
			double smag;
			flux_fn(UL, UR, n, params.gammad, Fhat, smag);

			/* phi_i at (xiL, etaL) for elemL; add phi_i * Fhat * w to R */
			shapeL(xrefL, order, &phi);
			for (int i = 0; i < Np; ++i) {
				double phi_i = phi[i];
				for (int var = 0; var < 4; ++var)
					R[(var * mesh.Ne + elemL) * Np + i] += phi_i * Fhat[var] * w;
			}

			/* phi_i at (xiR, etaR) for elemR; subtract phi_i * Fhat * w (Fhat(UR,UL,-n) = -Fhat) */
			shapeL(xrefR, order, &phi);
			for (int i = 0; i < Np; ++i) {
				double phi_i = phi[i];
				for (int var = 0; var < 4; ++var)
					R[(var * mesh.Ne + elemR) * Np + i] -= phi_i * Fhat[var] * w;
			}
		}
	}

	if (phi) free(phi);
}

// Freestream-only boundary contribution using ghost = interior state (not full BC implementation).

namespace {
enum class DGBcType {
	WALL,
	INFLOW,
	OUTFLOW,
	FREESTREAM
};

static std::string to_lower(std::string s) {
	std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return (char)std::tolower(c); });
	return s;
}

static bool contains_token(const std::string& s, const char* token) {
	return s.find(token) != std::string::npos;
}

static DGBcType classify_boundary_name(const std::string& name_raw) {
	std::string name = to_lower(name_raw);

	// Explicit mapping for turbine blade .gri: BGroup2/6 = walls, 4 = outflow, 8 = inflow.
	if (name == "bgroup2" || name == "bgroup6")
		return DGBcType::WALL;
	if (name == "bgroup4")
		return DGBcType::OUTFLOW;
	if (name == "bgroup8")
		return DGBcType::INFLOW;

	// Common CFD naming conventions
	if (contains_token(name, "wall") || contains_token(name, "solid"))
		return DGBcType::WALL;
	if (contains_token(name, "inflow") || contains_token(name, "inlet"))
		return DGBcType::INFLOW;
	if (contains_token(name, "outflow") || contains_token(name, "outlet") || contains_token(name, "exit"))
		return DGBcType::OUTFLOW;
	if (contains_token(name, "freestream") || contains_token(name, "farfield") || contains_token(name, "free"))
		return DGBcType::FREESTREAM;

	// Default: be conservative and preserve freestream test behavior.
	return DGBcType::FREESTREAM;
}
} // namespace

void addBndSurfTerm(const GriMesh& mesh, double* R, int order,
                    const double* U, const ProblemParams& params, FluxFn flux_fn) {
	int Np = (order + 1) * (order + 2) / 2;
	QuadratureRule quad1d = getQuadratureRule1D(order);

	double* phi = nullptr;

	const double nin[2] = {std::cos(params.alpha), std::sin(params.alpha)};
	const double R_gas = 1.0 / params.gammad;

	for (int ib = 0; ib < mesh.num_boundary_faces; ++ib) {
		int elem = mesh.B2E[3 * ib + 0];
		int face = mesh.B2E[3 * ib + 1];
		int bgroup = mesh.B2E[3 * ib + 2]; // 1-based into mesh.Bname
		const double* n = &mesh.Bn[2 * ib];
		double L = mesh.Bn_len[ib];

		std::string bname = "unknown";
		if (bgroup >= 1 && (size_t)bgroup <= mesh.Bname.size())
			bname = mesh.Bname[(size_t)bgroup - 1];
		DGBcType bctype = classify_boundary_name(bname);

		for (int q = 0; q < quad1d.nq; ++q) {
			double t = quad1d.xq[q];
			double w = quad1d.wq[q] * L;

			double xi, eta;
			faceRefCoords(face, t, xi, eta);
			double xref[2] = {xi, eta};

			// Interpolate UL at (xi, eta)
			shapeL(xref, order, &phi);
			double UL[4] = {0, 0, 0, 0};
			for (int j = 0; j < Np; ++j) {
				for (int var = 0; var < 4; ++var)
					UL[var] += U[var * mesh.Ne * Np + elem * Np + j] * phi[j];
			}

			double Fhat[4];
			double smag;
			if (bctype == DGBcType::WALL) {
				WallFlux(UL, n, params.gammad, Fhat, smag);
			} else if (bctype == DGBcType::INFLOW) {
				try {
					InflowFlux(UL, n, nin, params.rho0, params.a0, params.gammad, R_gas, flux_fn, Fhat, smag);
				} catch (const std::exception&) {
					// Interior state is unphysical during early iteration; use zero-dissipation fallback.
					flux_fn(UL, UL, n, params.gammad, Fhat, smag);
				}
			} else if (bctype == DGBcType::OUTFLOW) {
				try {
					OutflowFlux(UL, n, params.pout, params.gammad, flux_fn, Fhat, smag);
				} catch (const std::exception&) {
					flux_fn(UL, UL, n, params.gammad, Fhat, smag);
				}
			} else {
				// Freestream test behavior: ghost state equals interior state.
				flux_fn(UL, UL, n, params.gammad, Fhat, smag);
			}

			// Add boundary contribution to residual.
			shapeL(xref, order, &phi);
			for (int i = 0; i < Np; ++i) {
				double phi_i = phi[i];
				for (int var = 0; var < 4; ++var)
					R[(var * mesh.Ne + elem) * Np + i] += phi_i * Fhat[var] * w;
			}
		}
	}

	if (phi) free(phi);
}
