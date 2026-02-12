cmake_minimum_required(VERSION 3.10)
project(ysk_lib LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Recommended warnings (non-fatal)
if (MSVC)
  add_compile_options(/W4)
else()
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# Library sources
set(SOURCES
  src/Trapezoid.cpp
  src/Simpson13.cpp
  src/GaussianQuad.cpp
  src/GaussianQuad2D.cpp
  src/MonteCarlo.cpp
  src/TestFunctions.cpp
)

add_library(ysk_integration ${SOURCES})
target_include_directories(ysk_integration PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# 1D executable
add_executable(main apps/main.cpp)
target_link_libraries(main PRIVATE ysk_integration)

# 2D executable (extension)
add_executable(main2D apps/main2D.cpp)
target_link_libraries(main2D PRIVATE ysk_integration)
