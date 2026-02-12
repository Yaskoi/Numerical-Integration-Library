// ============================================
// include/solvers/GaussianQuadrature.hpp
// ============================================
#ifndef GAUSSIANQUADRATURE_HPP
#define GAUSSIANQUADRATURE_HPP

#include "Solver.hpp"
#include <vector>
#include <utility>

class GaussianQuadrature : public Solver {
public:
    GaussianQuadrature() = default;
    ~GaussianQuadrature() override = default;

    double integrate(const Function& f, double a, double b, int n) const override;

private:
    std::vector<std::pair<double, double>> getGaussPointsAndWeights(int n) const;
};

#endif // GAUSSIANQUADRATURE_HPP
