#ifndef TESTFUNCTIONS_HPP
#define TESTFUNCTIONS_HPP

#include "Function.hpp"
#include <cmath>

// x^2 * cos(x)
class QuadraticCosine : public Function {
public:
    double operator()(double x) const override {
        return x * x * std::cos(x);
    }
};

// x^10
class PowerTen : public Function {
public:
    double operator()(double x) const override {
        return std::pow(x, 10);
    }
};

// x^(-1/2)
class InverseSqrt : public Function {
public:
    double operator()(double x) const override {
        return 1.0 / std::sqrt(x);
    }
};

// log(x)
class Logarithm : public Function {
public:
    double operator()(double x) const override {
        return std::log(x);
    }
};

// 1/(x²+x^3+2)
class InvertedSumFunction : public Function {
public:
    double operator()(double x) const override {
        return 1.0/(x*x + x*x*x + 2.0);
    }
};

#endif