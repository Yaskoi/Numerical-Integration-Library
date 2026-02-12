// src/GaussianQuadrature.cpp

#include "solvers/GaussianQuad.hpp"
#include <stdexcept>
#include <cmath>

std::vector<std::pair<double, double>>
GaussianQuadrature::getGaussPointsAndWeights(int n) const {

    std::vector<std::pair<double, double>> points_weights;

    switch(n) {
        case 1:
            // n=1: Exact pour polynômes de degré ≤ 1
            points_weights.push_back({0.0, 2.0});
            break;

        case 2:
            // n=2: Exact pour polynômes de degré ≤ 3
            points_weights.push_back({-1.0/std::sqrt(3.0), 1.0});
            points_weights.push_back({ 1.0/std::sqrt(3.0), 1.0});
            break;

        case 3:
            // n=3: Exact pour polynômes de degré ≤ 5
            points_weights.push_back({-std::sqrt(3.0/5.0), 5.0/9.0});
            points_weights.push_back({ 0.0,                8.0/9.0});
            points_weights.push_back({ std::sqrt(3.0/5.0), 5.0/9.0});
            break;

        case 4:
            // n=4: Exact pour polynômes de degré ≤ 7
            {
                double x1 = std::sqrt(3.0/7.0 - (2.0/7.0)*std::sqrt(6.0/5.0));
                double x2 = std::sqrt(3.0/7.0 + (2.0/7.0)*std::sqrt(6.0/5.0));
                double w1 = (18.0 + std::sqrt(30.0)) / 36.0;
                double w2 = (18.0 - std::sqrt(30.0)) / 36.0;

                points_weights.push_back({-x2, w2});
                points_weights.push_back({-x1, w1});
                points_weights.push_back({ x1, w1});
                points_weights.push_back({ x2, w2});
            }
            break;

        case 5:
            // n=5: Exact pour polynômes de degré ≤ 9
            {
                double x1 = (1.0/3.0) * std::sqrt(5.0 - 2.0*std::sqrt(10.0/7.0));
                double x2 = (1.0/3.0) * std::sqrt(5.0 + 2.0*std::sqrt(10.0/7.0));
                double w0 = 128.0/225.0;
                double w1 = (322.0 + 13.0*std::sqrt(70.0)) / 900.0;
                double w2 = (322.0 - 13.0*std::sqrt(70.0)) / 900.0;

                points_weights.push_back({-x2, w2});
                points_weights.push_back({-x1, w1});
                points_weights.push_back({ 0.0, w0});
                points_weights.push_back({ x1, w1});
                points_weights.push_back({ x2, w2});
            }
            break;

        default:
            throw std::invalid_argument("Gaussian Quadrature: n doit être entre 1 et 5");
    }

    return points_weights;
}

double GaussianQuadrature::integrate(const Function& f, double a, double b, int n) const {
    // Vérification des paramètres
    if (n < 1) {
        throw std::invalid_argument("Gaussian Quadrature: n doit être ≥ 1");
    }

    if (a >= b) {
        throw std::invalid_argument("La borne inférieure doit être strictement inférieure à la borne supérieure");
    }

    // Si n ≤ 5: utilise directement n points de Gauss sur tout l'intervalle
    if (n <= 5) {
        // Récupérer les points et poids pour [-1, 1]
        auto points_weights = getGaussPointsAndWeights(n);

        // Changement de variable de [-1, 1] vers [a, b]
        double mid = (b + a) / 2.0;
        double half_length = (b - a) / 2.0;

        double sum = 0.0;
        for (const auto& [point, weight] : points_weights) {
            double x = mid + half_length * point;
            sum += weight * f(x);
        }

        return half_length * sum;
    }

    // Si n > 5: Gauss composite - divise en n sous-intervalle et on applique Gauss-3 points sur chaque sous-intervalle
    else {
        double h = (b - a) / n;  // Largeur de chaque sous-intervalle
        double total = 0.0;

        // Récupérer les points et poids pour Gauss-3
        auto points_weights = getGaussPointsAndWeights(3);

        // Pour chaque sous-intervalle [xi, xi+1]
        for (int i = 0; i < n; ++i) {
            double ai = a + i * h;      // Début du sous-intervalle
            double bi = ai + h;          // Fin du sous-intervalle

            double mid = (bi + ai) / 2.0;
            double half_length = (bi - ai) / 2.0;

            double sum = 0.0;
            for (const auto& [point, weight] : points_weights) {
                double x = mid + half_length * point;
                sum += weight * f(x);
            }

            total += half_length * sum;
        }

        return total;
    }
}