// src/Trapezoid.cpp
#include "solvers/Trapezoid.hpp"
#include "solvers/Trapezoid2D.hpp"
#include <stdexcept>

double Trapezoid::integrate(const Function& f, double a, double b, int n) const {
    // Vérification des paramètres
    if (n <= 0) {
        throw std::invalid_argument("Le nombre de sous-intervalles doit être positif");
    }

    if (a >= b) {
        throw std::invalid_argument("La borne inférieure doit être strictement inférieure à la borne supérieure");
    }

    // Calcul de la largeur de chaque sous-intervalle
    double h = (b - a) / n;

    // Initialisation avec les valeurs aux bornes
    // Formule: I ≈ h/2 * [f(a) + 2*Σf(xi) + f(b)]
    double sum = f(a) + f(b);

    // Ajout des termes intermédiaires (multipliés par 2)
    for (int i = 1; i < n; ++i) {
        double xi = a + i * h;
        sum += 2.0 * f(xi);
    }

    // Multiplication par h/2
    return (h / 2.0) * sum;
}

// IIMPLEMENTATION TRAPEZOID 2D

double Trapezoid2D::integrate(const Function2D& f,
                              double a, double b,
                              double c, double d,
                              int nx, int ny) const
{
    // Vérification des paramètres
    if (nx <= 0 || ny <= 0) {
        throw std::invalid_argument("Les nombres de subdivisions doivent être positifs");
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

    // Formule des trapèzes 2D:
    // ∫∫ f(x,y) dy dx ≈ (hx·hy/4) × [
    //   f(coins) + 2×f(bords_x) + 2×f(bords_y) + 4×f(intérieur)
    // ]

    double sum = 0.0;

    // Parcourir tous les points de la grille
    for (int i = 0; i <= nx; ++i) {
        double x = a + i * hx;

        for (int j = 0; j <= ny; ++j) {
            double y = c + j * hy;

            // Déterminer le coefficient multiplicateur
            double weight = 1.0;

            // Points intérieurs (ni sur le bord en x, ni sur le bord en y)
            if (i > 0 && i < nx && j > 0 && j < ny) {
                weight = 4.0;
            }
            // Points sur les bords (mais pas les coins)
            else if ((i == 0 || i == nx) && (j > 0 && j < ny)) {
                // Bord gauche ou droit, pas les coins
                weight = 2.0;
            }
            else if ((j == 0 || j == ny) && (i > 0 && i < nx)) {
                // Bord bas ou haut, pas les coins
                weight = 2.0;
            }
            // Sinon c'est un coin: weight = 1.0 (déjà initialisé)

            sum += weight * f(x, y);
        }
    }

    // Multiplication par (hx·hy/4)
    return (hx * hy / 4.0) * sum;
}