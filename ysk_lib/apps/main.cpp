// apps/main.cpp
#include <iostream>
#include <iomanip>
#include <cmath>
#include <functional>
#include <vector>
#include "Function.hpp"
#include "Solver.hpp"
#include "solvers/Trapezoid.hpp"
#include "solvers/Simpson13.hpp"
#include "solvers/GaussianQuad.hpp"
#include "solvers/MonteCarlo.hpp"
#include "functions/TestFunctions.hpp"

// Structure pour stocker les résultats
struct TestResult {
    std::string method_name;
    int n;
    double result;
    double exact;
    double error;
};

void printHeader(const std::string& title) {
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << title << "\n";
    std::cout << std::string(60, '=') << "\n\n";
}

void testFunction(const Function& f, double a, double b, double exact,
                  const std::string& test_name, const std::string& formula,
                  const std::vector<Solver*>& solvers,
                  const std::vector<int>& n_values,
                  const std::function<double(double)>& postprocess = nullptr) {

    printHeader(test_name + ": " + formula);
    std::cout << "Valeur exacte: " << exact << "\n\n";

    for (const auto* solver : solvers) {
        std::string solver_name = "Méthode inconnue";
        if (dynamic_cast<const Trapezoid*>(solver)) {
            solver_name = "Trapezoid";
        } else if (dynamic_cast<const Simpson13*>(solver)) {
            solver_name = "Simpson's 1/3 rule";
        } else if (dynamic_cast<const GaussianQuadrature*>(solver)) {
            solver_name = "Gaussian Quadrature";
        } else if (dynamic_cast<const MonteCarlo*>(solver)) {
            solver_name = "Monte-Carlo simulations";
        }

        std::cout << solver_name << ":\n";

        for (int n : n_values) {
            try {
                const double raw = solver->integrate(f, a, b, n);
                const double result = postprocess ? postprocess(raw) : raw;
                const double error = std::abs(result - exact);

                std::cout << "  n = " << std::setw(5) << n
                          << " : " << std::setw(15) << result
                          << " (erreur: " << std::scientific << error
                          << std::fixed << ")\n";
            } catch (const std::exception& e) {
                std::cout << "  n = " << std::setw(5) << n
                          << " : ERREUR - " << e.what() << "\n";
            }
        }
        std::cout << "\n";
    }
}

void compareConvergence(const Function& f, double a, double b, double exact,
                       const std::vector<Solver*>& solvers) {

    printHeader("COMPARAISON DE LA CONVERGENCE");

    std::vector<int> n_values = {2, 4, 10, 50, 100, 200, 500, 1000};

    std::cout << std::setw(8) << "n";
    for (const auto* solver : solvers) {
        std::string name = "Inconnu";
        if (dynamic_cast<const Trapezoid*>(solver)) name = "Trapezoid";
        else if (dynamic_cast<const Simpson13*>(solver)) name = "Simpson's 1/3";
        else if (dynamic_cast<const GaussianQuadrature*>(solver)) name = "Gaussian Quad";
        else if (dynamic_cast<const MonteCarlo*>(solver)) name = "Monte-Carlo";
        std::cout << std::setw(20) << name;
    }
    std::cout << "\n" << std::string(100, '-') << "\n";

    for (int n : n_values) {
        std::cout << std::setw(8) << n;
        for (const auto* solver : solvers) {
            try {
                double result = solver->integrate(f, a, b, n);
                double error = std::abs(result - exact);
                std::cout << std::setw(20) << std::scientific << error << std::fixed;
            } catch (...) {
                std::cout << std::setw(20) << "N/A";
            }
        }
        std::cout << "\n";
    }


    std::cout << "\nObservation: \n";
    std::cout << "  - Gaussian Quadrature: Le plus precis pour fonctions lisses\n";
    std::cout << "  - Simpson's 1/3 rule O(1/n^4): Tres bon compromis\n";
    std::cout << "  - Trapezoid O(1/n^2): Simple mais moins efficace\n";
    std::cout << "  - Monte-Carlo simulation O(1/n^(1/2)): Lent mais fonctionne en haute dimension\n";
}

int main() {
    std::cout << std::fixed << std::setprecision(10);

    printHeader("TESTS D'INTEGRATION NUMERIQUE");

    // Création des solveurs
    Trapezoid trapezoid;
    Simpson13 simpson13;
    GaussianQuadrature gauss;
    MonteCarlo montecarlo(42);  // Graine fixe pour reproductibilité

    std::vector<Solver*> solvers = {&trapezoid, &simpson13, &gauss, &montecarlo};

    // Valeurs de n à tester
    // Pour Gauss: n ≤ 5 → n points, n > 5 → n sous-intervalles avec Gauss-3
    // Pour Simpson: n doit être pair
    // Pour Trapèze: n quelconque
    std::vector<int> n_values = {2, 4, 10, 50, 100, 200, 500, 1000};

    // ========== TEST 1 : x² * cos(x) sur [0,1] ==========
    QuadraticCosine f1;
    double exact1 = 2*std::cos(1) - std::sin(1);
    testFunction(f1, 0.0, 1.0, exact1,
                 "Test 1", "QuadraticCosinuse",
                 solvers, n_values);

    // ========== TEST 2 : x^10 sur [0,1] ==========
    PowerTen f2;
    double exact2 = 1.0/11.0;
    testFunction(f2, 0.0, 1.0, exact2,
                 "Test 2", "PowerTen",
                 solvers, n_values);

    // ========== TEST 3 : x^(-1/2) sur [0,1] ==========
    InverseSqrt f3;
    double epsilon = 1e-6;
    // Solution 2 (découpage analytique):
    // ∫_0^1 f(x)dx = ∫_0^ε f(x)dx + ∫_ε^1 f(x)dx, avec ∫_0^ε x^{-1/2} dx = 2*sqrt(ε)
    const auto corr_inv_sqrt = [epsilon](double numeric_eps_1) {
        return numeric_eps_1 + 2.0 * std::sqrt(epsilon);
    };
    double exact3 = 2.0;  // Valeur exacte sur [0,1]
    testFunction(f3, epsilon, 1.0, exact3,
                 "Test 3", "InvertedSqrt",
                 solvers, n_values,
                 corr_inv_sqrt);

    // ========== TEST 4 : log(x) sur [0,1] ==========
    Logarithm f4;
    const auto corr_log = [epsilon](double numeric_eps_1) {
        return numeric_eps_1 + epsilon * (std::log(epsilon) - 1.0);
    };
    double exact4 = -1.0;  // Valeur exacte sur [0,1]
    testFunction(f4, epsilon, 1.0, exact4,
                 "Test 4", "Logarithm",
                 solvers, n_values,
                 corr_log);

    // ========== TEST 5 : 1/(x²+x³+2) sur [0,1] ==========
    InvertedSumFunction f5;
    double exact5 = 0.4047218685;
    testFunction(f5, 0.0, 1.0, exact5,
                 "Test 5", "Inverted Sum Function",
                 solvers, n_values);

    // ========== COMPARAISON DE CONVERGENCE ==========
    compareConvergence(f1, 0.0, 1.0, exact1, solvers);

    return 0;
}