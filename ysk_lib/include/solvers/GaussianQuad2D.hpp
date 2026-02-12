#ifndef GAUSSIANQUADRATURE2D_HPP
#define GAUSSIANQUADRATURE2D_HPP

#include "Solver2D.hpp"
#include <vector>
#include <utility>

class GaussianQuadrature2D : public Solver2D {
public:
    GaussianQuadrature2D() = default;
    ~GaussianQuadrature2D() override = default;

    double integrate(const Function2D& f,
                     double a, double b,
                     double c, double d,
                     int nx, int ny) const override;

private:
    // Helper pour récupérer les points/poids 1D (réutilisé pour le produit tensoriel)
    std::vector<std::pair<double, double>> getGaussPointsAndWeights(int n) const;
};

#endif // GAUSSIANQUADRATURE2D_HPP