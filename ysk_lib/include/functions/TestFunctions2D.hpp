#ifndef TESTFUNCTIONS2D_HPP
#define TESTFUNCTIONS2D_HPP

#include "Function2D.hpp"
#include <cmath>

class SquaredSum2D : public Function2D {
public:
    double operator()(double x, double y) const override {
        return x*x + y*y;
    }
};

class TrigProduct2D : public Function2D {
public:
    double operator()(double x, double y) const override {
        return std::sin(x) * std::cos(y);
    }
};

// Sur [0,2]×[0,2]
class Gaussian2D : public Function2D {
public:
    double operator()(double x, double y) const override {
        return std::exp(-(x*x + y*y));
    }
};

// f(x,y) = (x³ - 3x²y + xy² + y³)
class ComplexPolynomial2D : public Function2D {
public:
    double operator()(double x, double y) const override {
        return x*x*x - 3*x*x*y + x*y*y + y*y*y;
    }
};

class ComplexTrigo2D: public Function2D {
public:
    double operator()(double x, double y) const override {
        return x*x*std::cos(y) + y*y*std::sin(x);
    }
};

#endif //TESTFUNCTIONS2D_HPP
