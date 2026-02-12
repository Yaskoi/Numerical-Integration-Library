// include/Solver.hpp
#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "Function.hpp"

class Solver {
public:
    virtual ~Solver() = default;
    virtual double integrate(const Function& f, double a, double b, int n) const = 0;
};

#endif