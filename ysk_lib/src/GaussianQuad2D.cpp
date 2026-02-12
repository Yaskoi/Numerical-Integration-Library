#include "solvers/GaussianQuad2D.hpp"
#include <stdexcept>
#include <cmath>

/*Si n <=5 : On applique Gauss globalement ]
Si n > 5 : On découpe le domaine en $n_x \times n_y$ petits rectangles et on applique Gauss-3 (3 points) 
sur chaque petit rectangle (méthode composite).
*/

std::vector<std::pair<double, double>>
GaussianQuadrature2D::getGaussPointsAndWeights(int n) const {
    std::vector<std::pair<double, double>> points_weights;

    switch(n) {
        case 1:
            points_weights.push_back({0.0, 2.0});
            break;
        case 2:
            points_weights.push_back({-1.0/std::sqrt(3.0), 1.0});
            points_weights.push_back({ 1.0/std::sqrt(3.0), 1.0});
            break;
        case 3:
            points_weights.push_back({-std::sqrt(3.0/5.0), 5.0/9.0});
            points_weights.push_back({ 0.0,                8.0/9.0});
            points_weights.push_back({ std::sqrt(3.0/5.0), 5.0/9.0});
            break;
        case 4:
            {
                double x1 = std::sqrt(3.0/7.0 - (2.0/7.0)*std::sqrt(6.0/5.0));
                double x2 = std::sqrt(3.0/7.0 + (2.0/7.0)*std::sqrt(6.0/5.0));
                double w1 = (18.0 + std::sqrt(30.0)) / 36.0;
                double w2 = (18.0 - std::sqrt(30.0)) / 36.0;
                points_weights.push_back({-x2, w2}); points_weights.push_back({-x1, w1});
                points_weights.push_back({ x1, w1}); points_weights.push_back({ x2, w2});
            }
            break;
        case 5:
            {
                double x1 = (1.0/3.0) * std::sqrt(5.0 - 2.0*std::sqrt(10.0/7.0));
                double x2 = (1.0/3.0) * std::sqrt(5.0 + 2.0*std::sqrt(10.0/7.0));
                double w0 = 128.0/225.0;
                double w1 = (322.0 + 13.0*std::sqrt(70.0)) / 900.0;
                double w2 = (322.0 - 13.0*std::sqrt(70.0)) / 900.0;
                points_weights.push_back({-x2, w2}); points_weights.push_back({-x1, w1});
                points_weights.push_back({ 0.0, w0});
                points_weights.push_back({ x1, w1}); points_weights.push_back({ x2, w2});
            }
            break;
        default:
            throw std::invalid_argument("Gaussian Quadrature 2D: n doit être entre 1 et 5 pour le mode direct");
    }
    return points_weights;
}

double GaussianQuadrature2D::integrate(const Function2D& f,
                                       double a, double b,
                                       double c, double d,
                                       int nx, int ny) const {
    if (nx < 1 || ny < 1) throw std::invalid_argument("nx et ny doivent être >= 1");
    if (a >= b || c >= d) throw std::invalid_argument("Bornes invalides");

    // === CAS 1 : MODE DIRECT (nx, ny petits) ===
    // On utilise nx points en X et ny points en Y sur TOUT le domaine
    if (nx <= 5 && ny <= 5) {
        auto pw_x = getGaussPointsAndWeights(nx);
        auto pw_y = getGaussPointsAndWeights(ny);

        double mid_x = (b + a) / 2.0;
        double half_x = (b - a) / 2.0;
        double mid_y = (d + c) / 2.0;
        double half_y = (d - c) / 2.0;

        double sum = 0.0;

        // Double boucle sur les points de Gauss (Produit tensoriel)
        for (const auto& px : pw_x) {
            double u = px.first;  // point dans [-1, 1]
            double w_u = px.second;
            double x = mid_x + half_x * u;

            for (const auto& py : pw_y) {
                double v = py.first; // point dans [-1, 1]
                double w_v = py.second;
                double y = mid_y + half_y * v;

                sum += w_u * w_v * f(x, y);
            }
        }
        
        return half_x * half_y * sum;
    }

    
    else {
        double hx = (b - a) / nx;
        double hy = (d - c) / ny;
        double total_integral = 0.0;

        // On utilise Gauss-3 car c'est un bon compromis
        auto pw = getGaussPointsAndWeights(3);

        for (int i = 0; i < nx; ++i) {
            double xi = a + i * hx;
            double xi_next = xi + hx;
            double mid_x = (xi + xi_next) / 2.0;
            double half_x = hx / 2.0;

            for (int j = 0; j < ny; ++j) {
                double yj = c + j * hy;
                double yj_next = yj + hy;
                double mid_y = (yj + yj_next) / 2.0;
                double half_y = hy / 2.0;

                double local_sum = 0.0;

                for (const auto& px : pw) {
                    double x = mid_x + half_x * px.first;
                    for (const auto& py : pw) {
                        double y = mid_y + half_y * py.first;
                        local_sum += px.second * py.second * f(x, y);
                    }
                }
                total_integral += half_x * half_y * local_sum;
            }
        }
        return total_integral;
    }
}