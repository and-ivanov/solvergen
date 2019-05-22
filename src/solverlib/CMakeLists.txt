cmake_minimum_required(VERSION 3.12)

project(SolverLib)

FetchContent_Declare(vectorclass
    URL http://www.agner.org/optimize/vectorclass.zip
)

FetchContent_GetProperties(vectorclass)
if(NOT vectorclass_POPULATED)
    FetchContent_Populate(vectorclass)
endif()

add_library(solver_lib solverhelper.cpp solverhelper.h)

target_include_directories(
    solver_lib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> # hardcoded for build
    $<INSTALL_INTERFACE:.>
)
target_include_directories(solver_lib PUBLIC ${vectorclass_SOURCE_DIR})

