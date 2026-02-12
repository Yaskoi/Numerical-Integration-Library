// ============================================
// include/solvers/Simpson13.hpp
// ============================================
#ifndef SIMPSON13_HPP
#define SIMPSON13_HPP

#include "Solver.hpp"

class Simpson13 : public Solver {
public:

    Simpson13() = default;
    ~Simpson13() override = default;
    double integrate(const Function& f, double a, double b, int n) const override;
};

#endif // SIMPSON13_HPP