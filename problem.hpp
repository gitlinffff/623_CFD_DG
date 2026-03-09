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
 * Initialize all cells to uniform state from inlet conditions at given Mach number.
 * Uses isentropic relations: p/p0 = (1 + (g-1)/2*M^2)^(-g/(g-1)), etc.
 * Flow direction: u = M*a*cos(alpha), v = M*a*sin(alpha).
 *
 * U: output, size Ne*4, conserved state [rho, rho*u, rho*v, E]
 * Ne: number of elements
 * M: Mach number (e.g. 0.1 for steady init)
 */
void initialize_uniform(double* U, int Ne, int order, double M, const ProblemParams& params) {

#endif
