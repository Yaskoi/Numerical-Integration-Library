#ifndef FUNCTION2D_HPP
#define FUNCTION2D_HPP

class Function2D {
public:
    virtual ~Function2D() = default;
    virtual double operator()(double x, double y) const = 0;
};

#endif // FUNCTION2D_HPP