// include/solvers/Trapezoid.hpp
#ifndef TRAPEZOID_HPP
#define TRAPEZOID_HPP

#include "Solver.hpp"

class Trapezoid : public Solver {
public:

    Trapezoid() = default;
    ~Trapezoid() override = default;

    double integrate(const Function& f, double a, double b, int n) const override;
};

#endif // TRAPEZOID_HPP
