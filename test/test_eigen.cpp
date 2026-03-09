#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>

int main() {
    // 1. Test Dense Matrix (like your Jacobian or Physical Flux)
    Eigen::Matrix2d A;
    A << 1, 2, 
         3, 4;
    Eigen::Vector2d b(2, 1);
    Eigen::Vector2d x = A * b;

    std::cout << "--- Dense Test ---" << std::endl;
    std::cout << "Matrix A:\n" << A << "\nResult A*b:\n" << x << std::endl;

    // 2. Test Sparse Matrix (like your Global Mass Matrix)
    Eigen::SparseMatrix<double> S(2, 2);
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.push_back(Eigen::Triplet<double>(0, 0, 10.0));
    triplets.push_back(Eigen::Triplet<double>(1, 1, 20.0));
    S.setFromTriplets(triplets.begin(), triplets.end());

    std::cout << "\n--- Sparse Test ---" << std::endl;
    std::cout << "Sparse S (0,0): " << S.coeff(0,0) << std::endl;
    
    if (x(0) == 4 && x(1) == 10) {
        std::cout << "\nSUCCESS: Eigen logic is correct!" << std::endl;
    } else {
        std::cout << "\nFAILURE: Eigen logic error." << std::endl;
    }

    return 0;
}
