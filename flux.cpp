#include "flux.hpp"
#include "physics.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace {
const double eps = 1e-14;
const double kappa = 0.1;
}

void fluxRusanov(const double UL[4], const double UR[4], const double n[2], double gammad,
             double Fhat[4], double& smag) {
    double rhoL, uL, vL, pL, cL;
    double rhoR, uR, vR, pR, cR;
    consToPrim(UL, gammad, rhoL, uL, vL, pL, cL);
    consToPrim(UR, gammad, rhoR, uR, vR, pR, cR);

    double nx = n[0];
    double ny = n[1];

    double unL = uL * nx + vL * ny;
    double unR = uR * nx + vR * ny;

    double alpha = std::max(std::abs(unL) + cL, std::abs(unR) + cR);
    smag = alpha;

    double FL[4], FR[4];
    physicalFlux(UL, n, gammad, FL);
    physicalFlux(UR, n, gammad, FR);

    for (int k = 0; k < 4; ++k)
        Fhat[k] = 0.5 * (FL[k] + FR[k]) - 0.5 * alpha * (UR[k] - UL[k]);
}


void fluxHLLC(const double UL[4], const double UR[4], const double n[2], double gammad,
          double Fhat[4], double& smag) {
    double rhoL, uL, vL, pL, cL;
    double rhoR, uR, vR, pR, cR;
    consToPrim(UL, gammad, rhoL, uL, vL, pL, cL);
    consToPrim(UR, gammad, rhoR, uR, vR, pR, cR);

    double FL[4], FR[4];
    physicalFlux(UL, n, gammad, FL);
    physicalFlux(UR, n, gammad, FR);

    double nx = n[0];
    double ny = n[1];

    double unL = uL * nx + vL * ny;
    double utL = -uL * ny + vL * nx;
    double unR = uR * nx + vR * ny;
    double utR = -uR * ny + vR * nx;

    double SL = std::min(unL - cL, unR - cR);
    double SR = std::max(unL + cL, unR + cR);
    smag = std::max(std::abs(SL), std::abs(SR));

    double denom = rhoL * (SL - unL) - rhoR * (SR - unR);
    if (std::abs(denom) < 1e-14) denom = 1e-14;
    double s_star = (pR - pL + rhoL * unL * (SL - unL) - rhoR * unR * (SR - unR)) / denom;
    double p_star = pL + rhoL * (SL - unL) * (s_star - unL);

    if (SL >= 0) {
        for (int k = 0; k < 4; ++k)
            Fhat[k] = FL[k];
        return;
    }
    if (SR <= 0) {
        for (int k = 0; k < 4; ++k)
            Fhat[k] = FR[k];
        return;
    }

    if (s_star >= 0) {
        double sl = SL;
        double un = unL;
        double r_star = rhoL * (sl - un) / (sl - s_star);
        double E_L = UL[3];
        double E_star = ((sl - un) * E_L - pL * un + p_star * s_star) / (sl - s_star);
        double u_star = s_star * nx - utL * ny;
        double v_star = s_star * ny + utL * nx;

        Fhat[0] = FL[0] + sl * (r_star - UL[0]);
        Fhat[1] = FL[1] + sl * (r_star * u_star - UL[1]);
        Fhat[2] = FL[2] + sl * (r_star * v_star - UL[2]);
        Fhat[3] = FL[3] + sl * (E_star - UL[3]);
    } else {
        double sr = SR;
        double un = unR;
        double r_star = rhoR * (sr - un) / (sr - s_star);
        double E_R = UR[3];
        double E_star = ((sr - un) * E_R - pR * un + p_star * s_star) / (sr - s_star);
        double u_star = s_star * nx - utR * ny;
        double v_star = s_star * ny + utR * nx;
    
        Fhat[0] = FR[0] + sr * (r_star - UR[0]);
        Fhat[1] = FR[1] + sr * (r_star * u_star - UR[1]);
        Fhat[2] = FR[2] + sr * (r_star * v_star - UR[2]);
        Fhat[3] = FR[3] + sr * (E_star - UR[3]);
    }
}


void fluxROE(const double UL[4], const double UR[4], const double n[2], double gammad,
         double Fhat[4], double& smag) {
    double nx = n[0];
    double ny = n[1];

    double rhoL, uL, vL, pL, cL;
    double rhoR, uR, vR, pR, cR;

    consToPrim(UL, gammad, rhoL, uL, vL, pL, cL);
    consToPrim(UR, gammad, rhoR, uR, vR, pR, cR);

    double unL = uL * nx + vL * ny;
    double utL = -uL * ny + vL * nx;
    double unR = uR * nx + vR * ny;
    double utR = -uR * ny + vR * nx;

    double HL = (UL[3] + pL) / rhoL;
    double HR = (UR[3] + pR) / rhoR;

    double sqrt_rhoL = std::sqrt(rhoL);
    double sqrt_rhoR = std::sqrt(rhoR);
    double denom = sqrt_rhoL + sqrt_rhoR;
    if (denom < eps) denom = eps;

    double un_tilde = (sqrt_rhoL * unL + sqrt_rhoR * unR) / denom;
    double ut_tilde = (sqrt_rhoL * utL + sqrt_rhoR * utR) / denom;
    double H_tilde = (sqrt_rhoL * HL + sqrt_rhoR * HR) / denom;

    double c_tilde_sq = (gammad - 1.0) * (H_tilde - 0.5 * (un_tilde * un_tilde + ut_tilde * ut_tilde));
    if (c_tilde_sq < eps) c_tilde_sq = eps;
    double c_tilde = std::sqrt(c_tilde_sq);

    double lam1 = un_tilde - c_tilde;
    double lam2 = std::abs(un_tilde);
    double lam3 = std::abs(un_tilde);
    double lam4 = un_tilde + c_tilde;

    double delta = kappa * c_tilde;
    double abs_lam1 = std::abs(lam1);
    lam1 = (abs_lam1 < delta) ? (lam1 * lam1 + delta * delta) / (2.0 * delta) : abs_lam1;
    double abs_lam4 = std::abs(lam4);
    lam4 = (abs_lam4 < delta) ? (lam4 * lam4 + delta * delta) / (2.0 * delta) : abs_lam4;

    smag = std::max(std::abs(un_tilde - c_tilde), std::abs(un_tilde + c_tilde));

    double drho = rhoR - rhoL;
    double dun = unR - unL;
    double dp = pR - pL;
    double rho_tilde = std::sqrt(rhoL * rhoR);
    double c_tilde2 = c_tilde * c_tilde;

    double alpha1 = (dp - rho_tilde * c_tilde * dun) / (2.0 * c_tilde2);
    double alpha2 = drho - dp / c_tilde2;
    double alpha3 = (UR[2] * nx - UR[1] * ny) - (UL[2] * nx - UL[1] * ny) - ut_tilde * drho;
    double alpha4 = (dp + rho_tilde * c_tilde * dun) / (2.0 * c_tilde2);

    double l1 = lam1, l2 = lam2, l3 = lam3, l4 = lam4;
    double a1 = alpha1, a2 = alpha2, a3 = alpha3, a4 = alpha4;
    double un = un_tilde, ut = ut_tilde, ct = c_tilde, Ht = H_tilde;

    double d0 = l1 * a1 * 1.0 + l2 * a2 * 1.0 + l3 * a3 * 0.0 + l4 * a4 * 1.0;
    double d_n = l1 * a1 * (un - ct) + l2 * a2 * un + l3 * a3 * 0.0 + l4 * a4 * (un + ct);
    double d_t = l1 * a1 * ut + l2 * a2 * ut + l3 * a3 * 1.0 + l4 * a4 * ut;
    double d3 = l1 * a1 * (Ht - un * ct) + l2 * a2 * (0.5 * (un * un + ut * ut)) + l3 * a3 * ut + l4 * a4 * (Ht + un * ct);

    double D[4];
    D[0] = d0;
    D[1] = d_n * nx - d_t * ny;
    D[2] = d_n * ny + d_t * nx;
    D[3] = d3;

    double FL[4], FR[4];
    physicalFlux(UL, n, gammad, FL);
    physicalFlux(UR, n, gammad, FR);

    for (int k = 0; k < 4; ++k)
        Fhat[k] = 0.5 * (FL[k] + FR[k]) - 0.5 * D[k];
}


void fluxHLLE(const double UL[4], const double UR[4], const double n[2], double gammad,
          double Fhat[4], double& smag) {
    double rhoL, uL, vL, pL, cL;
    double rhoR, uR, vR, pR, cR;
    consToPrim(UL, gammad, rhoL, uL, vL, pL, cL);
    consToPrim(UR, gammad, rhoR, uR, vR, pR, cR);

    double nx = n[0];
    double ny = n[1];

    double unL = uL * nx + vL * ny;
    double unR = uR * nx + vR * ny;

    double SL = std::min(unL - cL, unR - cR);
    double SR = std::max(unL + cL, unR + cR);
    smag = std::max(std::abs(SL), std::abs(SR));

    double FL[4], FR[4];
    physicalFlux(UL, n, gammad, FL);
    physicalFlux(UR, n, gammad, FR);

    if (SL >= 0.0) {
        for (int k = 0; k < 4; ++k)
            Fhat[k] = FL[k];
        return;
    }
    if (SR <= 0.0) {
        for (int k = 0; k < 4; ++k)
            Fhat[k] = FR[k];
        return;
    }

    double dS = SR - SL;
    if (std::abs(dS) < eps) dS = eps;
    for (int k = 0; k < 4; ++k)
        Fhat[k] = (SR * FL[k] - SL * FR[k] + SL * SR * (UR[k] - UL[k])) / dS;
}


void WallFlux(const double UL[4], const double n[2], double gammad,
         double Fhat[4], double& smag){
    double rhoL, uL, vL, pL, cL;
    consToPrim(UL, gammad, rhoL, uL, vL, pL, cL);

    double ub[2];
    double pb;
    ub[0] = uL - (uL*n[0] + vL*n[1])*n[0];
    ub[1] = vL - (uL*n[0] + vL*n[1])*n[1];
    pb = (gammad - 1) * (UL[3] - 0.5*rhoL*(ub[0]*ub[0] + ub[1]*ub[1]));
    
    Fhat[0] = 0.0;
    Fhat[1] = pb*n[0];
    Fhat[2] = pb*n[1];
    Fhat[3] = 0.0;
    smag = std::abs(uL*n[0] + vL*n[1]) + cL;
}


void InflowFlux(const double UL[4], const double n[2], const double nin[2],
                double rho0, double a0, double gammad, double R,
                void (*FluxFunction)(const double*, const double*, const double*, double, double*, double&),
                double Fhat[4], double& smag){
    const double pt = rho0 * a0 * a0 / gammad;
    const double Tt = (a0 * a0) / (gammad * R);

    double rhoL, uL, vL, pL, cL;
    consToPrim(UL, gammad, rhoL, uL, vL, pL, cL);

    double Jplus = (uL*n[0] +vL*n[1]) + 2*cL/(gammad - 1);
    double dn = n[0]*nin[0] + n[1]*nin[1];
    
    double A = gammad*R*Tt*dn*dn - 0.5*(gammad - 1)*Jplus*Jplus;
    double B = 4*gammad*R*Tt*dn/(gammad - 1);
    double C = 4*gammad*R*Tt/((gammad-1)*(gammad-1)) - Jplus*Jplus;

    double delta = B*B - 4*A*C;
    if (delta < 0) {
        throw std::runtime_error("Inflow BC: negative discriminant");
    }

    double Mb2 = (-B + std::sqrt(delta))/(2*A);
    double Mb1 = (-B - std::sqrt(delta))/(2*A);
    double Mb;

    if (Mb1 >= 0 && Mb2 >= 0)
        Mb = std::min(Mb1, Mb2);
    else if (Mb1 >= 0)
            Mb = Mb1;
    else if (Mb2 >= 0)
            Mb = Mb2;
    else
            throw std::runtime_error("Inflow BC: no physical Mach root");
    
    const double Tb  = Tt / (1.0 + 0.5*(gammad - 1.0)*Mb*Mb);

    const double pb = pt * std::pow(Tb / Tt, gammad/(gammad - 1.0));

    const double rhob = pb / (R * Tb);
    const double cb   = std::sqrt(std::max(0.0, gammad*pb/rhob));

    const double ub = Mb * cb * nin[0];
    const double vb = Mb * cb * nin[1];

    double Ub[4];
    Ub[0] = rhob;
    Ub[1] = rhob * ub;
    Ub[2] = rhob * vb;
    Ub[3] = pb/(gammad - 1.0) + 0.5*rhob*(ub*ub + vb*vb);

    FluxFunction(UL, Ub, n, gammad, Fhat, smag);
}


// test use
void InflowFlux_compute_Ub(const double UL[4], const double n[2], const double nin[2],
                           double rho0, double a0, double gammad, double R, double Ub[4]) {
    const double pt = rho0 * a0 * a0 / gammad;
    const double Tt = (a0 * a0) / (gammad * R);
    double rhoL, uL, vL, pL, cL;
    consToPrim(UL, gammad, rhoL, uL, vL, pL, cL);
    double Jplus = (uL*n[0] + vL*n[1]) + 2*cL/(gammad - 1);
    double dn = n[0]*nin[0] + n[1]*nin[1];
    double A = gammad*R*Tt*dn*dn - 0.5*(gammad - 1)*Jplus*Jplus;
    double B = 4*gammad*R*Tt*dn/(gammad - 1);
    double C = 4*gammad*R*Tt/((gammad-1)*(gammad-1)) - Jplus*Jplus;
    double delta = B*B - 4*A*C;
    if (delta < 0) throw std::runtime_error("InflowFlux_compute_Ub: negative discriminant");
    double Mb1 = (-B - std::sqrt(delta))/(2*A);
    double Mb2 = (-B + std::sqrt(delta))/(2*A);
    double Mb = (Mb1 >= 0 && Mb2 >= 0) ? std::min(Mb1, Mb2) : (Mb1 >= 0 ? Mb1 : Mb2);
    if (Mb < 0) throw std::runtime_error("InflowFlux_compute_Ub: no physical Mach root");
    const double Tb = Tt / (1.0 + 0.5*(gammad - 1.0)*Mb*Mb);
    const double pb = pt * std::pow(Tb / Tt, gammad/(gammad - 1.0));
    const double rhob = pb / (R * Tb);
    const double cb = std::sqrt(std::max(0.0, gammad*pb/rhob));
    const double vb0 = Mb * cb * nin[0];
    const double vb1 = Mb * cb * nin[1];
    Ub[0] = rhob;
    Ub[1] = rhob * vb0;
    Ub[2] = rhob * vb1;
    Ub[3] = pb/(gammad - 1.0) + 0.5*rhob*(vb0*vb0 + vb1*vb1);
}

void OutflowFlux(const double UL[4], const double n[2],
                double pout, double gammad,
                void (*FluxFunction)(const double*, const double*, const double*, double, double*, double&),
                double Fhat[4], double& smag){

    double rhoL, uL, vL, pL, cL;
    consToPrim(UL, gammad, rhoL, uL, vL, pL, cL);

    if (rhoL <= 0)
        throw std::runtime_error("Outflow BC: negative density");
    
    double unL = uL*n[0] + vL*n[1];
    double Jplus = unL + 2*cL/(gammad - 1);

    double Splus = pL/std::pow(rhoL, gammad);
    if (Splus <= 0)
        throw std::runtime_error("Outflow BC: negative entropy");

    double pb = pout;
    double rhob = std::pow((pb/Splus),1.0/gammad);
    if (rhob <= 0)
        throw std::runtime_error("Outflow BC: negative boundary density");

    double cb = std::sqrt(gammad*pb/rhob);

    double unb = Jplus - 2*cb/(gammad - 1);
    double ub = uL - unL*n[0] + unb*n[0];
    double vb = vL - unL*n[1] + unb*n[1];


    double Ub[4];
    Ub[0] = rhob;
    Ub[1] = rhob * ub;
    Ub[2] = rhob * vb;
    Ub[3] = pb/(gammad - 1.0) + 0.5*rhob*(ub*ub + vb*vb);

    FluxFunction(UL, Ub, n, gammad, Fhat, smag);
}

// test use
void OutflowFlux_compute_Ub(const double UL[4], const double n[2], double pout,
                            double gammad, double Ub[4]) {
    double rhoL, uL, vL, pL, cL;
    consToPrim(UL, gammad, rhoL, uL, vL, pL, cL);
    double unL = uL*n[0] + vL*n[1];
    double Jplus = unL + 2*cL/(gammad - 1);
    double Splus = pL/std::pow(rhoL, gammad);
    double pb = pout;
    double rhob = std::pow((pb/Splus), 1.0/gammad);
    double cb = std::sqrt(gammad*pb/rhob);
    double unb = Jplus - 2*cb/(gammad - 1);
    double ub = uL - unL*n[0] + unb*n[0];
    double vb = vL - unL*n[1] + unb*n[1];
    Ub[0] = rhob;
    Ub[1] = rhob * ub;
    Ub[2] = rhob * vb;
    Ub[3] = pb/(gammad - 1.0) + 0.5*rhob*(ub*ub + vb*vb);
}
