#include "surf.hpp"
#include "mm.hpp"
#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstring>
#include <string>
#include <iostream>

extern "C" {
	void shapeL(double* xref, int p, double** pphi);
}

namespace {

void faceRefCoords(int face, double t, double& xi, double& eta) {
	switch (face) {
		case 0: xi = t;       eta = 0.0;     break;
		case 1: xi = 1.0-t;   eta = t;       break;
		case 2: xi = 0.0;     eta = 1.0-t;   break;
		default: xi = t;      eta = 0.0;     break;
	}
}

auto check_state = [](const char* tag,
                       int elem,
                       int face,
                       int qid,
                       const double Uq[4],
                       double gamma) {
    for (int m = 0; m < 4; ++m) {
        if (!std::isfinite(Uq[m])) {
            std::cerr << tag << " nonfinite U at elem=" << elem
                      << " face=" << face
                      << " q=" << qid
                      << " U=(" << Uq[0] << ", " << Uq[1] << ", "
                      << Uq[2] << ", " << Uq[3] << ")\n";
            throw std::runtime_error("nonfinite reconstructed state");
        }
    }

    double rho  = Uq[0];
    double rhou = Uq[1];
    double rhov = Uq[2];
    double E    = Uq[3];

    if (rho <= 0.0) {
        std::cerr << tag << " negative rho at elem=" << elem
                  << " face=" << face
                  << " q=" << qid
                  << " rho=" << rho
                  << " U=(" << Uq[0] << ", " << Uq[1] << ", "
                  << Uq[2] << ", " << Uq[3] << ")\n";
        throw std::runtime_error("negative reconstructed density");
    }

    double u = rhou / rho;
    double v = rhov / rho;
    double p = (gamma - 1.0) * (E - 0.5 * rho * (u * u + v * v));

    if (!std::isfinite(p) || p <= 0.0) {
        std::cerr << tag << " bad pressure at elem=" << elem
                  << " face=" << face
                  << " q=" << qid
                  << " p=" << p
                  << " rho=" << rho
                  << " E=" << E
                  << " u=" << u
                  << " v=" << v
                  << " U=(" << Uq[0] << ", " << Uq[1] << ", "
                  << Uq[2] << ", " << Uq[3] << ")\n";
        throw std::runtime_error("nonphysical reconstructed pressure");
    }
};

void facePhysCoords(const GriMesh& mesh, int k, int face, double t,
                   double& x, double& y) {
	int v0 = mesh.E[k * 3 + face];
	int v1 = mesh.E[k * 3 + (face + 1) % 3];
	double x0 = mesh.V[v0 * 2 + 0], y0 = mesh.V[v0 * 2 + 1];
	double x1 = mesh.V[v1 * 2 + 0], y1 = mesh.V[v1 * 2 + 1];
	x = x0 + t * (x1 - x0);
	y = y0 + t * (y1 - y0);
}

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

bool boundary_nodes_available(const GriMesh& mesh) {
	return (int)mesh.BedgeNodeOffset.size() == mesh.num_boundary_faces + 1;
}

void quadratic_edge_geom(const GriMesh& mesh, int n0, int n1, int nm, double t,
                         double& nx, double& ny, double& jac) {
	double x0 = mesh.V[n0 * 2 + 0], y0 = mesh.V[n0 * 2 + 1];
	double x1 = mesh.V[n1 * 2 + 0], y1 = mesh.V[n1 * 2 + 1];
	double xm = mesh.V[nm * 2 + 0], ym = mesh.V[nm * 2 + 1];

	double dN0 = 4.0 * t - 3.0;
	double dNm = 4.0 - 8.0 * t;
	double dN1 = 4.0 * t - 1.0;

	double dxdt = dN0 * x0 + dNm * xm + dN1 * x1;
	double dydt = dN0 * y0 + dNm * ym + dN1 * y1;
	jac = std::sqrt(dxdt * dxdt + dydt * dydt);
	if (jac < 1e-20) jac = 1e-20;

	double tx = dxdt / jac;
	double ty = dydt / jac;
	nx = ty;
	ny = -tx;
}

}

void addSurfTerm(const GriMesh& mesh, double* R, int order,
                const double* U, const ProblemParams& params, FluxFn flux_fn,
                std::vector<double>& sum_s) {
	int Np = (order + 1) * (order + 2) / 2;
	QuadratureRule quad1d = getQuadratureRule1D(order);

	#pragma omp parallel
	{
		double* phi = nullptr;

		#pragma omp for schedule(static)
		for (int iface = 0; iface < mesh.num_interior_faces; ++iface) {
			int elemL = mesh.I2E[4 * iface + 0];
			int faceL = mesh.I2E[4 * iface + 1];
			int elemR = mesh.I2E[4 * iface + 2];
			int faceR = mesh.I2E[4 * iface + 3];

			const double* n = &mesh.In[2 * iface];
			double max_smag_edge = 0.0;
			double L = mesh.In_len[iface];

			for (int q = 0; q < quad1d.nq; ++q) {
				double t = quad1d.xq[q];
				double w = quad1d.wq[q] * L;

				double xiL, etaL;
				faceRefCoords(faceL, t, xiL, etaL);

				double x_phys, y_phys;
				facePhysCoords(mesh, elemL, faceL, t, x_phys, y_phys);

				double xiR, etaR;
				physToRef(mesh, elemR, x_phys, y_phys, xiR, etaR);
				if (xiR < -0.1 || etaR < -0.1 || xiR + etaR > 1.1) {

					faceRefCoords(faceR, t, xiR, etaR);
				}

				double xrefL[2] = {xiL, etaL};
				shapeL(xrefL, order, &phi);
				double UL[4] = {0, 0, 0, 0};
				for (int j = 0; j < Np; ++j) {
					for (int var = 0; var < 4; ++var)
						UL[var] += U[var * mesh.Ne * Np + elemL * Np + j] * phi[j];
				}

				double xrefR[2] = {xiR, etaR};
				shapeL(xrefR, order, &phi);
				double UR[4] = {0, 0, 0, 0};
				for (int j = 0; j < Np; ++j) {
					for (int var = 0; var < 4; ++var)
						UR[var] += U[var * mesh.Ne * Np + elemR * Np + j] * phi[j];
				}

				check_state("UL", elemL, faceL, q, UL, params.gammad);
				check_state("UR", elemR, faceR, q, UR, params.gammad);

				double Fhat[4], smag_q;
				flux_fn(UL, UR, n, params.gammad, Fhat, smag_q);
				max_smag_edge = std::max(max_smag_edge, smag_q);

				shapeL(xrefL, order, &phi);
				for (int i = 0; i < Np; ++i) {
					double phi_i = phi[i];
					for (int var = 0; var < 4; ++var) {
						#pragma omp atomic
						R[(var * mesh.Ne + elemL) * Np + i] += phi_i * Fhat[var] * w;
					}
				}

				shapeL(xrefR, order, &phi);
				for (int i = 0; i < Np; ++i) {
					double phi_i = phi[i];
					for (int var = 0; var < 4; ++var) {
						#pragma omp atomic
						R[(var * mesh.Ne + elemR) * Np + i] -= phi_i * Fhat[var] * w;
					}
				}
			}
			#pragma omp atomic
			sum_s[elemL] += max_smag_edge * L;
			#pragma omp atomic
			sum_s[elemR] += max_smag_edge * L;
		}

		if (phi) free(phi);
	}
}

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

	if (name == "bgroup2" || name == "bgroup6")
		return DGBcType::WALL;
	if (name == "bgroup4")
		return DGBcType::OUTFLOW;
	if (name == "bgroup8")
		return DGBcType::INFLOW;

	if (contains_token(name, "wall") || contains_token(name, "solid"))
		return DGBcType::WALL;
	if (contains_token(name, "inflow") || contains_token(name, "inlet"))
		return DGBcType::INFLOW;
	if (contains_token(name, "outflow") || contains_token(name, "outlet") || contains_token(name, "exit"))
		return DGBcType::OUTFLOW;
	if (contains_token(name, "freestream") || contains_token(name, "farfield") || contains_token(name, "free"))
		return DGBcType::FREESTREAM;

	return DGBcType::FREESTREAM;
}
}

void addBndSurfTerm(const GriMesh& mesh, double* R, int order,
                    const double* U, const ProblemParams& params, FluxFn flux_fn,
											std::vector<double>& sum_s, bool in_ptb, const double t) {
	int Np = (order + 1) * (order + 2) / 2;
	QuadratureRule quad1d_linear = getQuadratureRule1D(order);
	int curved_quad_order = order + 1;
	if (curved_quad_order > 4) curved_quad_order = 4;
	QuadratureRule quad1d_curved = getQuadratureRule1D(curved_quad_order);

	const double nin[2] = {std::cos(params.alpha), std::sin(params.alpha)};
	const double R_gas = 1.0 / params.gammad;
	const bool has_bedge_nodes = boundary_nodes_available(mesh);

	#pragma omp parallel
	{
		double* phi = nullptr;

		#pragma omp for schedule(static)
			for (int ib = 0; ib < mesh.num_boundary_faces; ++ib) {
				int elem = mesh.B2E[3 * ib + 0];
				int face = mesh.B2E[3 * ib + 1];
				int bgroup = mesh.B2E[3 * ib + 2];
				double max_smag_edge = 0.0;
				double edge_measure = 0.0;

				int bedge_start = 0;
				int bedge_nnode = 0;
				if (has_bedge_nodes) {
					bedge_start = mesh.BedgeNodeOffset[(size_t)ib];
					bedge_nnode = mesh.BedgeNodeOffset[(size_t)ib + 1] - bedge_start;
				}
				bool is_quadratic_curve = has_bedge_nodes && (bedge_nnode == 3);
				const QuadratureRule& quad1d = is_quadratic_curve ? quad1d_curved : quad1d_linear;

				std::string bname = "unknown";
				if (bgroup >= 1 && (size_t)bgroup <= mesh.Bname.size())
					bname = mesh.Bname[(size_t)bgroup - 1];
				DGBcType bctype = classify_boundary_name(bname);

				for (int q = 0; q < quad1d.nq; ++q) {
					double tq = quad1d.xq[q];
					double n_q[2];
					double w = 0.0;
					if (is_quadratic_curve) {
						int n0 = mesh.BedgeNodes[(size_t)bedge_start + 0];
						int n1 = mesh.BedgeNodes[(size_t)bedge_start + 1];
						int nm = mesh.BedgeNodes[(size_t)bedge_start + 2];
						double jac;
						quadratic_edge_geom(mesh, n0, n1, nm, tq, n_q[0], n_q[1], jac);
						w = quad1d.wq[q] * jac;
					} else {
						n_q[0] = mesh.Bn[2 * ib + 0];
						n_q[1] = mesh.Bn[2 * ib + 1];
						w = quad1d.wq[q] * mesh.Bn_len[ib];
					}
					edge_measure += w;

					double xi, eta;
					faceRefCoords(face, tq, xi, eta);
					double xref[2] = {xi, eta};

				shapeL(xref, order, &phi);
				double UL[4] = {0, 0, 0, 0};
				for (int j = 0; j < Np; ++j) {
					for (int var = 0; var < 4; ++var)
						UL[var] += U[var * mesh.Ne * Np + elem * Np + j] * phi[j];
				}

				check_state("UL_bnd", elem, face, q, UL, params.gammad);

					double Fhat[4], smag_q;
					if (bctype == DGBcType::WALL) {
						WallFlux(UL, n_q, params.gammad, Fhat, smag_q);
					}
					else if (bctype == DGBcType::INFLOW) {
						double rho0_in = params.rho0;
						if (in_ptb) {
							double x_phys, y_phys;
							facePhysCoords(mesh, elem, face, tq, x_phys, y_phys);
							double y_rot = y_phys;
							double ystator = y_rot + params.Vrot * t;
							double eta = ystator/params.delta_y - std::floor(ystator/params.delta_y) - 0.5;
							double fac = 1.0 - params.fwake * std::exp(-eta * eta / (2.0 * params.delta_wake * params.delta_wake));
							rho0_in = params.rho0 * fac;
						}
						try {
							InflowFlux(UL, n_q, nin, rho0_in, params.a0, params.gammad, R_gas, flux_fn, Fhat, smag_q);
						} catch (const std::exception&) {
							flux_fn(UL, UL, n_q, params.gammad, Fhat, smag_q);
						}
					}
					else if (bctype == DGBcType::OUTFLOW) {
						try {
							OutflowFlux(UL, n_q, params.pout, params.gammad, flux_fn, Fhat, smag_q);
						} catch (const std::exception&) {
							flux_fn(UL, UL, n_q, params.gammad, Fhat, smag_q);
						}
					}
					else {
						flux_fn(UL, UL, n_q, params.gammad, Fhat, smag_q);
					}
					max_smag_edge = std::max(max_smag_edge, smag_q);

				shapeL(xref, order, &phi);
				for (int i = 0; i < Np; ++i) {
					double phi_i = phi[i];
					for (int var = 0; var < 4; ++var) {
						#pragma omp atomic
						R[(var * mesh.Ne + elem) * Np + i] += phi_i * Fhat[var] * w;
					}
				}
				}
				#pragma omp atomic
				sum_s[elem] += max_smag_edge * edge_measure;
			}

		if (phi) free(phi);
	}
}
