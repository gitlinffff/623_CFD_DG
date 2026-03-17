#include "physics.hpp"
#include <cmath>
#include <algorithm>

void consToPrim(const double U[4], double gammad,
                  double& rho, double& u, double& v, double& p, double& c) {
    const double rho_floor = 1e-12;
    const double p_floor = 1e-12;

    rho = std::max(U[0], rho_floor);
    double mom_x = U[1];
    double mom_y = U[2];
    double E = U[3];

    u = mom_x / rho;
    v = mom_y / rho;

    p = (gammad - 1.0) * (E - 0.5 * rho * (u * u + v * v));
    p = std::max(p, p_floor);

    c = std::sqrt(gammad * p / rho);
}

void primToCons(double rho, double u, double v, double p, double gammad, double U[4]) {
    double E = p / (gammad - 1.0) + 0.5 * rho * (u * u + v * v);
    U[0] = rho;
    U[1] = rho * u;
    U[2] = rho * v;
    U[3] = E;
}

double getp0(double rho0, double a0, double gammad) {
    return rho0 * a0 * a0 / gammad;
}

void isentropic_prim_from_M(double rho0, double p0, double gammad, double M, double alpha,
                            double& rho, double& u, double& v, double& p) {

    double fac = 1.0 + 0.5 * (gammad - 1.0) * M * M;
    p = p0 / std::pow(fac, gammad / (gammad - 1.0));
    rho = rho0 / std::pow(fac, 1.0 / (gammad - 1.0));
    double a = std::sqrt(gammad * p / rho);
    u = M * a * std::cos(alpha);
    v = M * a * std::sin(alpha);
}

double total_temperature(double rho, double u, double v, double p, double gammad){

    const double eps = 1e-14;
    double a2 = gammad * p / (rho + eps);
    a2 = std::max(a2, eps);

    double M2 = (u*u + v*v) / a2;
    double fac = 1.0 + 0.5*(gammad-1.0)*M2;

    return gammad * (p / (rho + eps)) * fac;
}

double total_pressure(double rho, double u, double v, double p, double gammad) {

    const double eps = 1e-14;
    double a2 = gammad * p / (rho + eps);
    if (a2 < eps) a2 = eps;
    double M2 = (u * u + v * v) / a2;
    double fac = 1.0 + 0.5 * (gammad - 1.0) * M2;
    return p * std::pow(fac, gammad / (gammad - 1.0));
}

void physicalFlux(const double U[4], const double n[2], double gammad, double F[4]) {
    double nx = n[0];
    double ny = n[1];

    double rho, u, v, p, c;
    consToPrim(U, gammad, rho, u, v, p, c);
    double E = U[3];

    double un = u * nx + v * ny;

    F[0] = rho * un;
    F[1] = rho * u * un + p * nx;
    F[2] = rho * v * un + p * ny;
    F[3] = (E + p) * un;
}

Eigen::Matrix<double, 4, 2> vec_phyFlux(const double U[4], double gammad) {
	double rho, u, v, p, c;
	consToPrim(U, gammad, rho, u, v, p, c);
	double E = U[3];

	Eigen::Matrix<double, 4, 2> F;

	F(0, 0) = rho * u;
	F(1, 0) = rho * u * u + p;
	F(2, 0) = rho * u * v;
	F(3, 0) = (E + p) * u;

	F(0, 1) = rho * v;
	F(1, 1) = rho * v * u;
	F(2, 1) = rho * v * v + p;
	F(3, 1) = (E + p) * v;

	return F;
}
