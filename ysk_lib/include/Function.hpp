// include/Function.hpp
#ifndef FUNCTION_HPP
#define FUNCTION_HPP

class Function {
public:
    virtual ~Function() = default;
    virtual double operator()(double x) const = 0;
};

#endif

//end Function.hpp