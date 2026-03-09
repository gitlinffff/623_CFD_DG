#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include "quad.hpp"

// Helper to check near-equality for doubles
bool is_near(double a, double b, double tol = 1e-12) {
    return std::abs(a - b) < tol;
}

void test_partition_of_unity(int order) {
    std::cout << "Testing Order p=" << order << "..." << std::endl;
    
    std::vector<BasisEval> phiq = computePhiQ(order);
    QuadratureRule quad = getQuadratureRule(order);
    int Np = (order + 1) * (order + 2) / 2;

    for (int q = 0; q < quad.nq; ++q) {
        double sum_phi = 0.0;
        double sum_dphi_dxi = 0.0;
        double sum_dphi_deta = 0.0;

        for (int i = 0; i < Np; ++i) {
            int idx = i * quad.nq + q;
            sum_phi += phiq[idx].phi;
            sum_dphi_dxi += phiq[idx].dphi_dxi;
            sum_dphi_deta += phiq[idx].dphi_deta;
        }

        // 1. Sum of basis functions must be 1.0
        if (!is_near(sum_phi, 1.0)) {
            std::cerr << "Fail: p=" << order << " quad point " << q 
                      << " sum(phi) = " << sum_phi << std::endl;
            exit(1);
        }

        // 2. Sum of gradients must be 0.0
        if (!is_near(sum_dphi_dxi, 0.0) || !is_near(sum_dphi_deta, 0.0)) {
            std::cerr << "Fail: p=" << order << " quad point " << q 
                      << " sum(grad) != 0" << std::endl;
            exit(1);
        }
    }
    std::cout << "  - Partition of Unity: OK" << std::endl;
}

void test_integral_of_constant(int order) {
    // Integrating 1.0 over the reference triangle (Area = 0.5)
    std::vector<BasisEval> phiq = computePhiQ(order);
    QuadratureRule quad = getQuadratureRule(order);
    int Np = (order + 1) * (order + 2) / 2;

    double area = 0.0;
    // For a nodal basis, the integral of 1.0 is the sum of integrals of phi_i
    for (int q = 0; q < quad.nq; ++q) {
        // Just sum the weights (since sum(phi) is 1)
        area += quad.wq[q];
    }

    if (!is_near(area, 0.5)) {
        std::cerr << "Fail: p=" << order << " Area = " << area << " (expected 0.5)" << std::endl;
        exit(1);
    }
    std::cout << "  - Reference Area Integration: OK" << std::endl;
}

int main() {
    try {
        for (int p = 0; p <= 3; ++p) {
            test_partition_of_unity(p);
            test_integral_of_constant(p);
        }
        std::cout << "\nALL QUADRATURE TESTS PASSED!" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
