#include "surf.hpp"
#include "geometry.hpp"
#include <algorithm>
#include <cmath>
#include <cctype>
#include <string>
#include <iostream>
#include <stdexcept>

extern "C" {
    void shapeL(double* xref, int p, double** pphi);
}

namespace {

void check_state(const char* tag,
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
}

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

    return DGBcType::FREESTREAM;
}

} // namespace

void addSurfTerm(const GriMesh& mesh, double* R, int order,
                 const double* U, const ProblemParams& params, FluxFn flux_fn,
                 std::vector<double>* sum_s) {
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

            double max_smag_edge = 0.0;
            double face_length = 0.0;

            for (int q = 0; q < quad1d.nq; ++q) {
                double t = quad1d.xq[q];

                double n[2];
                double x_phys, y_phys;
                double ds_dt = face_metric_normal_point(mesh, elemL, faceL, t, n, x_phys, y_phys);
                double w = quad1d.wq[q] * ds_dt;
                face_length += w;

                // Left reference point on local face.
                double xiL, etaL;
                faceRefCoords(faceL, t, xiL, etaL);
                double xrefL[2] = {xiL, etaL};

                // Right reference point from physical mapping; periodic fallback uses local face map.
                double xiR, etaR;
                bool ok = phys_to_ref(mesh, elemR, x_phys, y_phys, xiR, etaR);
                if (!ok || xiR < -0.1 || etaR < -0.1 || xiR + etaR > 1.1) {
                    faceRefCoords(faceR, t, xiR, etaR);
                }
                double xrefR[2] = {xiR, etaR};

                // Interpolate UL.
                shapeL(xrefL, order, &phi);
                double UL[4] = {0, 0, 0, 0};
                for (int j = 0; j < Np; ++j) {
                    for (int var = 0; var < 4; ++var)
                        UL[var] += U[var * mesh.Ne * Np + elemL * Np + j] * phi[j];
                }

                // Interpolate UR.
                shapeL(xrefR, order, &phi);
                double UR[4] = {0, 0, 0, 0};
                for (int j = 0; j < Np; ++j) {
                    for (int var = 0; var < 4; ++var)
                        UR[var] += U[var * mesh.Ne * Np + elemR * Np + j] * phi[j];
                }

                check_state("UL", elemL, faceL, q, UL, params.gammad);
                check_state("UR", elemR, faceR, q, UR, params.gammad);

                // Numerical flux Fhat(UL, UR, nL_outward).
                double Fhat[4], smag_q;
                flux_fn(UL, UR, n, params.gammad, Fhat, smag_q);
                max_smag_edge = std::max(max_smag_edge, smag_q);

                // Add to left residual.
                shapeL(xrefL, order, &phi);
                for (int i = 0; i < Np; ++i) {
                    double phi_i = phi[i];
                    for (int var = 0; var < 4; ++var) {
                        #pragma omp atomic
                        R[(var * mesh.Ne + elemL) * Np + i] += phi_i * Fhat[var] * w;
                    }
                }

                // Subtract from right residual for conservation.
                shapeL(xrefR, order, &phi);
                for (int i = 0; i < Np; ++i) {
                    double phi_i = phi[i];
                    for (int var = 0; var < 4; ++var) {
                        #pragma omp atomic
                        R[(var * mesh.Ne + elemR) * Np + i] -= phi_i * Fhat[var] * w;
                    }
                }
            }

            if (sum_s != nullptr) {
                #pragma omp atomic
                (*sum_s)[elemL] += max_smag_edge * face_length;
                #pragma omp atomic
                (*sum_s)[elemR] += max_smag_edge * face_length;
            }
        }

        if (phi) free(phi);
    }
}

void addBndSurfTerm(const GriMesh& mesh, double* R, int order,
                    const double* U, const ProblemParams& params, FluxFn flux_fn,
                    std::vector<double>* sum_s) {
    int Np = (order + 1) * (order + 2) / 2;
    QuadratureRule quad1d = getQuadratureRule1D(order);

    const double nin[2] = {std::cos(params.alpha), std::sin(params.alpha)};
    const double R_gas = 1.0 / params.gammad;

    #pragma omp parallel
    {
        double* phi = nullptr;

        #pragma omp for schedule(static)
        for (int ib = 0; ib < mesh.num_boundary_faces; ++ib) {
            int elem = mesh.B2E[3 * ib + 0];
            int face = mesh.B2E[3 * ib + 1];
            int bgroup = mesh.B2E[3 * ib + 2]; // 1-based index in mesh.Bname

            double max_smag_edge = 0.0;
            double face_length = 0.0;

            std::string bname = "unknown";
            if (bgroup >= 1 && (size_t)bgroup <= mesh.Bname.size())
                bname = mesh.Bname[(size_t)bgroup - 1];
            DGBcType bctype = classify_boundary_name(bname);

            for (int q = 0; q < quad1d.nq; ++q) {
                double t = quad1d.xq[q];
                double xi, eta;
                faceRefCoords(face, t, xi, eta);
                double xref[2] = {xi, eta};

                double n[2];
                double x_phys, y_phys;
                double ds_dt = face_metric_normal_point(mesh, elem, face, t, n, x_phys, y_phys);
                double w = quad1d.wq[q] * ds_dt;
                face_length += w;

                // Interpolate UL at boundary quadrature point.
                shapeL(xref, order, &phi);
                double UL[4] = {0, 0, 0, 0};
                for (int j = 0; j < Np; ++j) {
                    for (int var = 0; var < 4; ++var)
                        UL[var] += U[var * mesh.Ne * Np + elem * Np + j] * phi[j];
                }

                check_state("UL_bnd", elem, face, q, UL, params.gammad);

                double Fhat[4], smag_q;
                if (bctype == DGBcType::WALL) {
                    WallFlux(UL, n, params.gammad, Fhat, smag_q);
                } else if (bctype == DGBcType::INFLOW) {
                    try {
                        InflowFlux(UL, n, nin, params.rho0, params.a0, params.gammad, R_gas, flux_fn, Fhat, smag_q);
                    } catch (const std::exception&) {
                        flux_fn(UL, UL, n, params.gammad, Fhat, smag_q);
                    }
                } else if (bctype == DGBcType::OUTFLOW) {
                    try {
                        OutflowFlux(UL, n, params.pout, params.gammad, flux_fn, Fhat, smag_q);
                    } catch (const std::exception&) {
                        flux_fn(UL, UL, n, params.gammad, Fhat, smag_q);
                    }
                } else {
                    // Freestream-test behavior.
                    flux_fn(UL, UL, n, params.gammad, Fhat, smag_q);
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

            if (sum_s != nullptr) {
                #pragma omp atomic
                (*sum_s)[elem] += max_smag_edge * face_length;
            }
        }

        if (phi) free(phi);
    }
}
