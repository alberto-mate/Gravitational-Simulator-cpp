cmake_minimum_required(VERSION 3.10)
project(PROYECTORENDIMIENTO)

set(CMAKE_CXX_STANDARD 17)
add_compile_options(-Wall -Wextra -Wno-deprecated -pedantic -pedantic-errors -DCMAKE_BUILD_TYPE=Release -Ofast -march=native)

add_executable(sim-aos sim-aos.cpp)
add_executable(sim-soa sim-soa.cpp)
add_executable(sim-aos-opti sim-aos-opti.cpp)
add_executable(sim-soa-opti sim-soa-opti.cpp)
