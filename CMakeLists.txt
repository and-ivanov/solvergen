cmake_minimum_required(VERSION 3.12)

project(SolverGen)

include(FetchContent)

FetchContent_Declare(hrdata
    GIT_REPOSITORY https://github.com/am-ivanov/hrdata.git
    GIT_TAG 9bb7a8334065801f503bdbc4b11ce03921fc0bee
)

FetchContent_GetProperties(hrdata)
if(NOT hrdata_POPULATED)
    FetchContent_Populate(hrdata)
    add_subdirectory(${hrdata_SOURCE_DIR} ${hrdata_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

FetchContent_Declare(llvm
    URL https://releases.llvm.org/8.0.0/llvm-8.0.0.src.tar.xz
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/llvm
)

FetchContent_GetProperties(llvm)
if(NOT llvm_POPULATED)
    FetchContent_Populate(llvm)
endif()

FetchContent_Declare(clang
    URL https://releases.llvm.org/8.0.0/cfe-8.0.0.src.tar.xz
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/llvm/tools/clang
)

FetchContent_GetProperties(clang)
if(NOT clang_POPULATED)
    FetchContent_Populate(clang)
endif()


function(add_llvm)
	set(CMAKE_BUILD_TYPE Release)
	message("INSIDE:" ${CMAKE_BUILD_TYPE})
	set(LLVM_ENABLE_RTTI ON) # otherwise it is impossible to inherit from llvm classess
	set(BUILD_SHARED_LIBS ON) # enables shared shibraries to link with
	#set(LLVM_LINK_LLVM_DYLIB ON)
	#set(LLVM_BUILD_LLVM_DYLIB ON)
	add_subdirectory(${llvm_SOURCE_DIR} ${llvm_BINARY_DIR} EXCLUDE_FROM_ALL)
	#add_subdirectory(${llvm_SOURCE_DIR} ${llvm_BINARY_DIR})
endfunction()

message("BEFORE:" ${CMAKE_BUILD_TYPE})
add_llvm()
message("AFTER:" ${CMAKE_BUILD_TYPE})

# create compilation database (compile_commands.json) to process later
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_executable(solver_generator src/generator/generator.cpp)

# because LLVM itself doesn't define include dirs
target_include_directories(solver_generator PRIVATE
	${llvm_SOURCE_DIR}/include
	${clang_SOURCE_DIR}/include
	${llvm_BINARY_DIR}/tools/clang/include
	${llvm_BINARY_DIR}/include
)

# it makes headers visible for libTooling based tools in build dir
add_dependencies(solver_generator clang-headers)

target_link_libraries(solver_generator PRIVATE
    LLVMOption
    LLVMSupport
    clangAST
    clangASTMatchers
    clangBasic
    clangDriver
    clangFrontend
    clangLex
    clangRewrite
    clangRewriteFrontend
    clangSerialization
    clangStaticAnalyzerFrontend
    clangTooling
    clangSema
)

set_target_properties(solver_generator PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${llvm_BINARY_DIR}/bin"
)

add_subdirectory(src/solverlib)

function(add_solver solver)
    # add entry to compilation database
    add_executable(${solver}-dummy "${CMAKE_CURRENT_SOURCE_DIR}/src/solvers/${solver}.cpp")
    target_link_libraries(${solver}-dummy hrdata solver_lib)

    add_custom_command(
        OUTPUT ${solver}.gen.cpp

        COMMAND ${llvm_BINARY_DIR}/bin/solver_generator
            ${CMAKE_CURRENT_SOURCE_DIR}/src/solvers/${solver}.cpp
            -o ${solver}.gen.cpp
            -p ${CMAKE_CURRENT_BINARY_DIR}

        # next line forces regeneration if generator or source are changed
        DEPENDS solver_generator ${solver}-dummy

        VERBATIM
    )

    add_executable(${solver} ${solver}.gen.cpp)
    target_link_libraries(${solver} hrdata solver_lib)
    if (${CMAKE_BUILD_TYPE} EQUAL "Debug")
        target_compile_definitions(${solver} PRIVATE _GLIBCXX_DEBUG)
    endif()
endfunction()


add_solver(gc_elastic_2d)
add_solver(gc_elastic_2d_from_le2d)
add_solver(fd_scalar_wave_2d)



