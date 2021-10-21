cmake_minimum_required(VERSION 3.10)
project(PROYECTORENDIMIENTO)

set(CMAKE_CXX_STANDARD 17)
add_compile_options(-Wall -Wextra -Wno-deprecated -DCMAKE_BUILD_TYPE=Release -pedantic -pedantic-errors -Ofast -march=native)

add_executable(sim-aos sim-aos.cpp)
add_executable(sim-soa sim-soa.cpp)