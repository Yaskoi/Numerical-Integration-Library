#ifndef SOLVER2D_HPP
#define SOLVER2D_HPP

#include "Function2D.hpp"

class Solver2D {
public:
    virtual ~Solver2D() = default;

    virtual double integrate(const Function2D& f,
                            double a, double b,
                            double c, double d,
                            int nx, int ny) const = 0;
};

#endif // SOLVER2D_HPP
