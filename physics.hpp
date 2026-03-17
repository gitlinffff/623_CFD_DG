#ifndef PHYSICS_HPP
#define PHYSICS_HPP

#include <Eigen/Dense>

void consToPrim(const double U[4], double gammad,
                  double& rho, double& u, double& v, double& p, double& c);

void primToCons(double rho, double u, double v, double p, double gammad, double U[4]);

double getp0(double rho0, double a0, double gammad);

void isentropic_prim_from_M(double rho0, double p0, double gammad, double M, double alpha,
                            double& rho, double& u, double& v, double& p);

double total_temperature(double rho, double u, double v, double p, double gammad);

double total_pressure(double rho, double u, double v, double p, double gammad);

void physicalFlux(const double U[4], const double n[2], double gammad, double Fout[4]);

Eigen::Matrix<double, 4, 2> vec_phyFlux(const double U[4], double gammad);

#endif
