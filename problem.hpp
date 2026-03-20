#ifndef PROBLEM_HPP
#define PROBLEM_HPP

/** Problem parameters for turbine blade passage (proj.pdf). */
struct ProblemParams {
    double rho0;   /**< Inlet stagnation density (e.g. 1) */
    double a0;     /**< Inlet stagnation speed of sound (e.g. 1) */
    double gammad;  /**< Ratio of specific heats (1.4) */
    double alpha;  /**< Inlet angle of attack [rad] (50 deg) */
    double pout;   /**< Outflow static pressure (0.7 * p0) */

    /** Unsteady wake inflow (proj.pdf): rho0(eta)=rho0*[1-fwake*exp(-eta^2/(2*delta^2))] */
    double Vrot;       /**< Rotor speed = a0 */
    double delta_y;    /**< Stator pitch, 18mm = 0.018 if mesh in m */
    double fwake;     /**< Wake deficit strength, 0.1 */
    double delta_wake; /**< Wake Gaussian width, 0.1 */

    /** Default: rho0=1, a0=1, gammad=1.4, alpha=50deg, pout=0.7*p0 */
    ProblemParams();
};

/** Stagnation pressure p0 from params */
double getp0(const ProblemParams& p);

/**
 * Initialize all DG DOFs to a uniform freestream state based on the inlet
 * stagnation conditions and a fixed Mach number M=0.1 (as in Project 2).
 *
 * U layout (DG): U[var * Ne * Np + k * Np + i].
 */
void initialize_uniform(double* U, int Ne, int order, const ProblemParams& params);

#endif
