// src/MonteCarlo.cpp

#include "solvers/MonteCarlo.hpp"
#include "solvers/MonteCarlo2D.hpp"
#include <stdexcept>
#include <cmath>

MonteCarlo::MonteCarlo(unsigned int seed)
    : generator(seed == 0 ? std::random_device{}() : seed),
      distribution(0.0, 1.0) {
}

double MonteCarlo::integrate(const Function& f, double a, double b, int n) const {
    // Vérification des paramètres
    if (n <= 0) {
        throw std::invalid_argument("Le nombre de simulations doit être positif");
    }

    if (a >= b) {
        throw std::invalid_argument("La borne inférieure doit être strictement inférieure à la borne supérieure");
    }

    // Méthode Monte-Carlo simple:
    // 1. Générer n points aléatoires uniformément sur [a, b]
    // 2. Évaluer f en chaque point
    // 3. Calculer la moyenne: I ≈ (b-a) × moyenne(f(Xi))

    double sum = 0.0;
    double length = b - a;

    // Générer n points aléatoires et évaluer f
    for (int i = 0; i < n; ++i) {
        // Générer un point aléatoire uniforme sur [0,1]
        double u = distribution(generator);

        // Transformer en point sur [a, b]
        double x = a + length * u;

        // Accumuler f(x)
        sum += f(x);
    }

    // Calculer l'intégrale: (b-a) × moyenne
    return length * (sum / n);
}

double MonteCarlo::estimateError(const Function& f, double a, double b, int n, double result) const {
    if (n <= 1) {
        return std::numeric_limits<double>::infinity();
    }

    // Calculer la variance de f sur [a,b]
    double sum_squared = 0.0;
    double length = b - a;
    double mean = result / length;  // Moyenne de f

    for (int i = 0; i < n; ++i) {
        double u = distribution(generator);
        double x = a + length * u;
        double fx = f(x);
        sum_squared += (fx - mean) * (fx - mean);
    }

    double variance = sum_squared / (n - 1);

    // Erreur standard: σ/√n × (b-a)
    return length * std::sqrt(variance / n);
}

// IMPLÉMENTATION MONTE-CARLO 2D

MonteCarlo2D::MonteCarlo2D(unsigned int seed)
    : generator(seed == 0 ? std::random_device{}() : seed),
      distribution(0.0, 1.0) {
}

double MonteCarlo2D::integrate(const Function2D& f,
                               double a, double b,
                               double c, double d,
                               int nx, int ny) const {
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

    // Nombre total de points aléatoires
    int n = nx * ny;

    // Calcul de l'aire du rectangle
    double area = (b - a) * (d - c);

    // Méthode Monte-Carlo 2D:
    // 1. Générer n points (x,y) aléatoires uniformément sur [a,b]×[c,d]
    // 2. Évaluer f en chaque point
    // 3. Calculer: I ≈ Aire × moyenne(f(xi, yi))

    double sum = 0.0;

    for (int i = 0; i < n; ++i) {
        // Générer un point aléatoire (x, y) dans [a,b] × [c,d]
        double u = distribution(generator);  // u ∈ [0,1]
        double v = distribution(generator);  // v ∈ [0,1]

        // Transformer en coordonnées du domaine
        double x = a + (b - a) * u;
        double y = c + (d - c) * v;

        // Accumuler f(x, y)
        sum += f(x, y);
    }

    // Calculer l'intégrale: Aire × moyenne
    return area * (sum / n);
}