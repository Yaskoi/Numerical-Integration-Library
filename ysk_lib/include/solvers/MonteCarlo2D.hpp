#ifndef MONTECARLO2D_HPP
#define MONTECARLO2D_HPP

#include "Solver2D.hpp"
#include <random>

class MonteCarlo2D : public Solver2D {
public:
    explicit MonteCarlo2D(unsigned int seed = 0);
    ~MonteCarlo2D() override = default;

    double integrate(const Function2D& f,
                     double a, double b,
                     double c, double d,
                     int nx, int ny) const override;

private:
    mutable std::mt19937 generator;
    mutable std::uniform_real_distribution<double> distribution;
};

#endif //MONTECARLO2D_HPP
