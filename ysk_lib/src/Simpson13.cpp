// src/Simpson13.cpp

#include "solvers/Simpson13.hpp"
#include "solvers/Simpson132D.hpp"
#include <stdexcept>

double Simpson13::integrate(const Function& f, double a, double b, int n) const {
    // Vérification des paramètres
    if (n <= 0) {
        throw std::invalid_argument("Le nombre de sous-intervalles doit être positif");
    }

    if (n % 2 != 0) {
        throw std::invalid_argument("Simpson 1/3 nécessite un nombre PAIR de sous-intervalles");
    }

    if (a >= b) {
        throw std::invalid_argument("La borne inférieure doit être strictement inférieure à la borne supérieure");
    }

    // Calcul de la largeur de chaque sous-intervalle
    double h = (b - a) / n;

    // Formule de Simpson 1/3:
    // I ≈ h/3 * [f(a) + 4*Σf(x_impair) + 2*Σf(x_pair) + f(b)]

    // Initialisation avec les bornes
    double sum = f(a) + f(b);

    // Somme des termes aux indices impairs (multipliés par 4)
    // i = 1, 3, 5, ..., n-1
    for (int i = 1; i < n; i += 2) {
        double xi = a + i * h;
        sum += 4.0 * f(xi);
    }

    // Somme des termes aux indices pairs (multipliés par 2)
    // i = 2, 4, 6, ..., n-2
    for (int i = 2; i < n; i += 2) {
        double xi = a + i * h;
        sum += 2.0 * f(xi);
    }

    // Multiplication par h/3
    return (h / 3.0) * sum;
}

// IMPLEMENTATION SIMPSON'S 1/3 2D

double Simpson132D::getSimpsonWeight(int i, int n) const {
    // Poids de Simpson 1D:
    // - i = 0 ou i = n : poids = 1 (bords)
    // - i impair : poids = 4
    // - i pair (mais pas bords) : poids = 2

    if (i == 0 || i == n) {
        return 1.0;
    } else if (i % 2 == 1) {
        return 4.0;
    } else {
        return 2.0;
    }
}

double Simpson132D::integrate(const Function2D& f,
                           double a, double b,
                           double c, double d,
                           int nx, int ny) const {
    // Vérification des paramètres
    if (nx <= 0 || ny <= 0) {
        throw std::invalid_argument("Les nombres de subdivisions doivent être positifs");
    }

    if (nx % 2 != 0 || ny % 2 != 0) {
        throw std::invalid_argument("Simpson 2D nécessite nx ET ny PAIRS");
    }

    if (a >= b) {
        throw std::invalid_argument("La borne inférieure en x doit être < à la borne supérieure");
    }

    if (c >= d) {
        throw std::invalid_argument("La borne inférieure en y doit être < à la borne supérieure");
    }

    // Calcul des pas
    double hx = (b - a) / nx;
    double hy = (d - c) / ny;

    // Formule de Simpson 2D:
    // ∫∫ f(x,y) dy dx ≈ (hx·hy/9) × Σᵢ Σⱼ wᵢ·wⱼ·f(xᵢ, yⱼ)

    double sum = 0.0;

    // Parcourir tous les points de la grille
    for (int i = 0; i <= nx; ++i) {
        double x = a + i * hx;
        double wx = getSimpsonWeight(i, nx);

        for (int j = 0; j <= ny; ++j) {
            double y = c + j * hy;
            double wy = getSimpsonWeight(j, ny);

            // Le poids 2D est le produit des poids 1D
            double weight = wx * wy;

            sum += weight * f(x, y);
        }
    }

    // Multiplication par (hx·hy/9)
    return (hx * hy / 9.0) * sum;
}