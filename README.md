# ysk_lib — Numerical Integration Library

A C++ library for numerical integration of 1D and 2D functions, developed as part of a C++ training project.

## Methods Implemented

### 1D Integration
- **Trapezoidal Rule** — O(1/n²)
- **Simpson's 1/3 Rule** — O(1/n⁴)
- **Gaussian Quadrature** — Exact for polynomials of degree ≤ 2n-1
- **Monte-Carlo** — O(1/√n)

### 2D Integration (Extension)
- **Trapezoid 2D**
- **Simpson 2D**
- **Gaussian Quadrature 2D**
- **Monte-Carlo 2D**

## Project Structure
```
ysk_lib/
├── CMakeLists.txt
├── include/
│   ├── Function.hpp
│   ├── Solver.hpp
│   ├── Function2D.hpp
│   ├── Solver2D.hpp
│   ├── solvers/
│   │   ├── Trapezoid.hpp
│   │   ├── Simpson13.hpp
│   │   ├── GaussianQuadrature.hpp
│   │   ├── MonteCarlo.hpp
│   │   ├── Trapezoid2D.hpp
│   │   ├── Simpson2D.hpp
│   │   ├── GaussianQuadrature2D.hpp
│   │   └── MonteCarlo2D.hpp
│   └── functions/
│       ├── TestFunctions.hpp
│       └── TestFunctions2D.hpp
├── src/
│   ├── Trapezoid.cpp
│   ├── Simpson13.cpp
│   ├── GaussianQuadrature.cpp
│   ├── GaussianQuadrature2D.cpp
│   ├── MonteCarlo.cpp
│   ├── TestFunctions.cpp
└── apps/
    ├── main.cpp
    └── main2D.cpp
```

## Build
```bash
mkdir build && cd build
cmake ..
make
./main      # 1D tests
./main2D    # 2D tests
```

## License

This project is licensed under the MIT License.
