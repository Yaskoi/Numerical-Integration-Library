// apps/main2D.cpp

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include "Function2D.hpp"
#include "Solver2D.hpp"
#include "solvers/Trapezoid2D.hpp"
#include "solvers/Simpson132D.hpp"
#include "solvers/MonteCarlo2D.hpp"
#include "functions/TestFunctions2D.hpp"
#include "solvers/GaussianQuad2D.hpp"

void printHeader(const std::string& title) {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << title << "\n";
    std::cout << std::string(70, '=') << "\n\n";
}

void testFunction2D(const Function2D& f,
                    double a, double b, double c, double d,
                    double exact,
                    const std::string& test_name,
                    const std::string& formula,
                    const std::vector<Solver2D*>& solvers,
                    const std::vector<int>& n_values) {

    printHeader(test_name + ": " + formula);
    std::cout << "Domaine: [" << a << "," << b << "] × ["
              << c << "," << d << "]\n";
    std::cout << "Valeur exacte: " << exact << "\n\n";

    for (const auto* solver : solvers) {
        std::string solver_name = "Methode inconnue";
        if (dynamic_cast<const Trapezoid2D*>(solver)) {
            solver_name = "Trapezoid 2D";
        } else if (dynamic_cast<const Simpson132D*>(solver)) {
            solver_name = "Simpson13 2D";
        } else if (dynamic_cast<const GaussianQuadrature2D*>(solver)) {
            solver_name = "Gaussian Quad 2D";
        } else if (dynamic_cast<const MonteCarlo2D*>(solver)) {
            solver_name = "Monte-Carlo 2D";
        }

        std::cout << solver_name << ":\n";

        for (int n : n_values) {
            try {
                double result = solver->integrate(f, a, b, c, d, n, n);
                double error = std::abs(result - exact);

                std::cout << "  n = " << std::setw(4) << n
                          << " (total: " << std::setw(7) << (n+1)*(n+1) << " pts)"
                          << " : " << std::setw(15) << result
                          << " (erreur: " << std::scientific << error
                          << std::fixed << ")\n";
            } catch (const std::exception& e) {
                std::cout << "  n = " << std::setw(4) << n
                          << " : ERREUR - " << e.what() << "\n";
            }
        }
        std::cout << "\n";
    }
}

void compareConvergence(const Function2D& f,
                       double a, double b, double c, double d,
                       double exact,
                       const std::vector<Solver2D*>& solvers) {

    printHeader("COMPARAISON DE LA CONVERGENCE 2D");

    std::vector<int> n_values = {10, 50, 100, 200, 500, 1000};
    std::vector<int> n_simulations = {500, 1000, 2000, 5000, 7500, 10000};

    std::cout << std::setw(6) << "n" << std::setw(10) << "Nb pts";
    for (const auto* solver : solvers) {
        std::string name = "Inconnu";
        if (dynamic_cast<const Trapezoid2D*>(solver)) name = "Trapezoid";
        else if (dynamic_cast<const Simpson132D*>(solver)) name = "Simpson's 1/3";
        else if (dynamic_cast<const GaussianQuadrature2D*>(solver)) name = "Gaussian";
        else if (dynamic_cast<const MonteCarlo2D*>(solver)) name = "MonteCarlo";
        std::cout << std::setw(18) << name;
    }
    std::cout << "\n" << std::string(80, '-') << "\n";

    for (int n : n_values) {
        int total_pts = (n + 1) * (n + 1);
        std::cout << std::setw(6) << n << std::setw(10) << total_pts;

        for (const auto* solver : solvers) {
            try {
                double result = solver->integrate(f, a, b, c, d, n, n);
                double error = std::abs(result - exact);
                std::cout << std::setw(18) << std::scientific << error << std::fixed;
            } catch (...) {
                std::cout << std::setw(18) << "N/A";
            }
        }
        std::cout << "\n";
    }

    std::cout << "\nObservation:\n";
    std::cout << "  - Simpson converge le plus vite: O(1/n^4)\n";
    std::cout << "  - Trapezoid: O(1/n^2)\n";
    std::cout << "  - Monte-Carlo: O(1/n)\n";
}

int main() {
    std::cout << std::fixed << std::setprecision(12);

    printHeader("TESTS D'INTEGRATION NUMERIQUE 2D - FONCTIONS AVANCEES");

    // Création des solveurs
    Trapezoid2D trapezoid;
    Simpson132D simpson;
    GaussianQuadrature2D gauss;
    MonteCarlo2D montecarlo(42);  // Graine fixe pour reproductibilité

    std::vector<Solver2D*> solvers = {&trapezoid, &simpson, &gauss, &montecarlo};

    // Valeurs de n à tester (n pair pour Simpson)
    std::vector<int> n_values = {4, 10, 50, 100, 200, 500, 1000};

    // TEST 1 : x² + y² sur [0,1]×[0,1]
    SquaredSum2D f1;
    testFunction2D(f1, 0.0, 1.0, 0.0, 1.0, 2.0/3.0,
                   "Test 1", "x^2 + y^2",
                   solvers, n_values);

    // TEST 2 : sin(x)*cos(y) sur [0,π]×[0,π]
    TrigProduct2D f2;
    testFunction2D(f2, 0.0, M_PI, 0.0, M_PI, 0.0,
                   "Test 2", "sin(x)*cos(y)",
                   solvers, n_values);

    // TEST 3 : e^(-(x²+y²)) sur [0,2]×[0,2]
    Gaussian2D f3;
    Simpson132D simpson_ref;
    double exact3 = simpson_ref.integrate(f3, 0.0, 2.0, 0.0, 2.0, 100, 100);
    testFunction2D(f3, 0.0, 2.0, 0.0, 2.0, exact3,
                   "Test 3", "f(x,y) = e^(-(x^2+y^2)) [Gaussienne 2D]",
                   solvers, n_values);

    // TEST 4 : Polynôme complexe sur [0,1]×[0,1]
    ComplexPolynomial2D f4;
    testFunction2D(f4, 0.0, 1.0, 0.0, 1.0, 1.0/6.0,
                   "Test 4", "f(x,y) = x^3 - 3x^2y + xy^2 + y^3",
                   solvers, n_values);

    // TEST 5 : Trigonométrie complexe sur [0,π]×[0,π]
    ComplexTrigo2D f5;
    double exact5 = 2.0 * std::pow(M_PI, 3) / 3.0;
    testFunction2D(f5, 0.0, M_PI, 0.0, M_PI, exact5,
                   "Test 5", "f(x,y) = (x^2)cos(y) + (y^2)sin(x)",
                   solvers, n_values);

    // ANALYSE DE CONVERGENCE sur fonction polynomiale
    compareConvergence(f4, 0.0, 1.0, 0.0, 1.0, 1.0/6.0, solvers);

    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "TESTS 2D TERMINES - " << solvers.size() << " methodes testees\n";
    std::cout << std::string(70, '=') << "\n";

    return 0;
}