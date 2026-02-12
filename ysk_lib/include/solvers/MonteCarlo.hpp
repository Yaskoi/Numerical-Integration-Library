// ============================================
// include/solvers/MonteCarlo.hpp
// ============================================
#ifndef MONTECARLO_HPP
#define MONTECARLO_HPP

#include "Solver.hpp"
#include <random>

class MonteCarlo : public Solver {
public:
    explicit MonteCarlo(unsigned int seed = 0);

    ~MonteCarlo() override = default;
    double integrate(const Function& f, double a, double b, int n) const override;

    double estimateError(const Function& f, double a, double b, int n, double result) const;

private:
    mutable std::mt19937 generator;  // Générateur Mersenne Twister
    mutable std::uniform_real_distribution<double> distribution;
};

#endif // MONTECARLO_HPP