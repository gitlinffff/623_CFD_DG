#ifndef PROBLEM_HPP
#define PROBLEM_HPP

struct ProblemParams {
    double rho0;
    double a0;
    double gammad;
    double alpha;
    double pout;

    double Vrot;
    double delta_y;
    double fwake;
    double delta_wake;

    ProblemParams();
};

double getp0(const ProblemParams& p);

void initialize_uniform(double* U, int Ne, int order, const ProblemParams& params);

#endif
