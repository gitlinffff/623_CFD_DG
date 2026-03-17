#ifndef FLUX_HPP
#define FLUX_HPP

void fluxRusanov(const double UL[4], const double UR[4], const double n[2], double gammad,
             double Fhat[4], double& smag);

void fluxHLLC(const double UL[4], const double UR[4], const double n[2], double gammad,
          double Fhat[4], double& smag);

void fluxROE(const double UL[4], const double UR[4], const double n[2], double gammad,
         double Fhat[4], double& smag);

void fluxHLLE(const double UL[4], const double UR[4], const double n[2], double gammad,
          double Fhat[4], double& smag);

void WallFlux(const double UL[4], const double n[2], double gammad,
              double Fhat[4], double& smag);

void InflowFlux(const double UL[4], const double n[2], const double nin[2],
                double rho0, double a0, double gammad, double R,
                void (*FluxFunction)(const double*, const double*, const double*, double, double*, double&),
                double Fhat[4], double& smag);

void OutflowFlux(const double UL[4], const double n[2], double pout, double gammad,
                 void (*FluxFunction)(const double*, const double*, const double*, double, double*, double&),
                 double Fhat[4], double& smag);

void InflowFlux_compute_Ub(const double UL[4], const double n[2], const double nin[2],
                           double rho0, double a0, double gammad, double R, double Ub[4]);

void OutflowFlux_compute_Ub(const double UL[4], const double n[2], double pout,
                            double gammad, double Ub[4]);

#endif
