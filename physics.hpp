#ifndef PHYSICS_HPP
#define PHYSICS_HPP

/** U = [rho, rho*u, rho*v, E], p = (g-1)*(E - rho*(u^2+v^2)/2), c = sqrt(g*p/rho) */
void consToPrim(const double U[4], double gammad,
                  double& rho, double& u, double& v, double& p, double& c);

/** E = p/(g-1) + rho*(u^2+v^2)/2 */
void primToCons(double rho, double u, double v, double p, double gammad, double U[4]);

/** Stagnation pressure: p0 = rho0 * a0^2 / gammad */
double getp0(double rho0, double a0, double gammad);

/**
 * Isentropic primitive state from stagnation (rho0,p0), Mach M, flow angle alpha [rad].
 * p/p0 = (1 + (g-1)/2 * M^2)^(-g/(g-1)),  rho/rho0 = (1 + (g-1)/2 * M^2)^(-1/(g-1))
 * u = M*a*cos(alpha), v = M*a*sin(alpha), a = sqrt(gammad*p/rho)
 */
void isentropic_prim_from_M(double rho0, double p0, double gammad, double M, double alpha,
                            double& rho, double& u, double& v, double& p);

/** Total temperature: Tt/T = 1 + (g-1)/2 * M^2, M^2 = (u^2+v^2)/a^2 */
double total_temperature(double rho, double u, double v, double p, double gammad);

/** Total pressure: pt/p = (1 + (g-1)/2 * M^2)^(g/(g-1)) */
double total_pressure(double rho, double u, double v, double p, double gammad);

/** F·n: [rho*un, rho*u*un+p*nx, rho*v*un+p*ny, (E+p)*un], un = u·n */
void physicalFlux(const double U[4], const double n[2], double gammad, double Fout[4]);

Eigen::Matrix<double, 4, 2> vec_phyFlux(const double U[4], double gammad)

#endif
