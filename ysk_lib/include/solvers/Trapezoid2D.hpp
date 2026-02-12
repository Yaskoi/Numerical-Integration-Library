#ifndef TRAPEZOID2D_HPP
#define TRAPEZOID2D_HPP

#include "Solver2D.hpp"

class Trapezoid2D: public Solver2D {
public:
    Trapezoid2D() = default;
    ~Trapezoid2D() override = default;

    double integrate(const Function2D& f,
                     double a, double b,
                     double c, double d,
                     int nx, int ny) const override;
};

#endif //TRAPEZOID2D_HPP
