cmake_minimum_required(VERSION 3.10)
project(davidfolly)
enable_language(Fortran)
set(CMAKE_C_STANDARD 99)
set(SOURCE main.c get_red_mass.c globals.c read_write.c get_integral.c bg_solver.c get_pot.c get_pot.h bg_solver.h get_ho.c get_ho.h bg_solver_2.c bg_solver_2.h)

add_executable(${PROJECT_NAME} ${SOURCE})
target_link_libraries(${PROJECT_NAME} gmosh)
add_library(gmosh gmosh.f90)