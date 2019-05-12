cmake_minimum_required(VERSION 3.12)

project(SolverLib)

add_library(solver_lib solverhelper.cpp solverhelper.h)
target_include_directories(
    solver_lib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> # hardcoded for build
    $<INSTALL_INTERFACE:.>
)
