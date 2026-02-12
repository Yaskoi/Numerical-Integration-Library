#ifndef SIMPSON132D_HPP
#define SIMPSON132D_HPP

#include "Solver2D.hpp"

class Simpson132D : public Solver2D {
public:
    Simpson132D() = default;
    ~Simpson132D() override = default;

    double integrate(const Function2D& f,
                     double a, double b,
                     double c, double d,
                     int nx, int ny) const override;

private:
    double getSimpsonWeight(int i, int n) const;
};

#endif //SIMPSON132D_HPP
