##############################################
#              Project Setup                 #
##############################################
cmake_minimum_required(VERSION 3.24.0)
project(NnPDE_Study LANGUAGES CXX C)

include(./cmake-common.cmake)
include_directories(include)

set(dir ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${dir}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${dir}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${dir}/bin")

find_package(fmt)

##############################################
#              My Libraries                  #
##############################################
add_library(prm_handler
        src/EllipticParamsHandler.cpp
        src/NavierStokesParamsHandler.cpp
        src/ParabolicParamsHandler.cpp
        src/DNParamHandler.cpp
        include/PDEParamsHandler.hpp
)
deal_ii_setup_target(prm_handler)
target_link_libraries(prm_handler fmt::fmt)

add_library(elliptic_lib src/EllipticSolver.cpp include/EllipticSolver.hpp)
deal_ii_setup_target(elliptic_lib)
target_link_libraries(elliptic_lib fmt::fmt prm_handler)

add_library(s_elliptic_lib src/SimpleGridEllipticSolver.cpp include/SimpleGridEllipticSolver.hpp)
deal_ii_setup_target(s_elliptic_lib)
target_link_libraries(s_elliptic_lib fmt::fmt prm_handler)

add_library(parabolic_lib src/ParabolicSolver.cpp include/ParabolicSolver.hpp)
deal_ii_setup_target(parabolic_lib)
target_link_libraries(parabolic_lib fmt::fmt prm_handler)

add_library(s_parabolic_lib src/SimpleGridParabolicSolver.cpp include/SimpleGridParabolicSolver.hpp)
deal_ii_setup_target(s_parabolic_lib)
target_link_libraries(s_parabolic_lib fmt::fmt prm_handler)

add_library(stokes_lib src/StokesSolver.cpp include/StokesSolver.hpp include/NavierStokesHelpers.hpp)
deal_ii_setup_target(stokes_lib)
target_link_libraries(stokes_lib fmt::fmt prm_handler)

add_library(unsteady_stokes_lib src/UnsteadyStokesSolver.cpp include/UnsteadyStokesSolver.hpp include/NavierStokesHelpers.hpp)
deal_ii_setup_target(unsteady_stokes_lib)
target_link_libraries(unsteady_stokes_lib fmt::fmt prm_handler)

add_library(dn_lib src/DNSolver.cpp src/DNEllipticSolver.cpp include/DNSolver.hpp include/DNEllipticSolver.hpp)
deal_ii_setup_target(dn_lib)
target_link_libraries(dn_lib fmt::fmt prm_handler)

##############################################
#                 My Tests                   #
##############################################
add_executable(test_ell test/test_elliptic.cpp)
deal_ii_setup_target(test_ell)
target_link_libraries(test_ell fmt::fmt prm_handler elliptic_lib )

add_executable(test_par test/test_parabolic.cpp)
deal_ii_setup_target(test_par)
target_link_libraries(test_par fmt::fmt prm_handler parabolic_lib)

add_executable(test_stokes test/test_stokes.cpp)
deal_ii_setup_target(test_stokes)
target_link_libraries(test_stokes fmt::fmt prm_handler stokes_lib)

add_executable(test_unsteady test/test_unsteady.cpp)
deal_ii_setup_target(test_unsteady)
target_link_libraries(test_unsteady fmt::fmt prm_handler unsteady_stokes_lib)

add_executable(test_dn test/test_dn.cpp)
deal_ii_setup_target(test_dn)
target_link_libraries(test_dn fmt::fmt prm_handler dn_lib)

add_executable(test_s_ell test/test_simple_elliptic.cpp)
deal_ii_setup_target(test_s_ell)
target_link_libraries(test_s_ell fmt::fmt prm_handler s_elliptic_lib)

add_executable(test_s_par test/test_simple_parabolic.cpp)
deal_ii_setup_target(test_s_par)
target_link_libraries(test_s_par fmt::fmt prm_handler s_parabolic_lib)

add_executable(test_par_conv test/conv_test_parabolic.cpp)
deal_ii_setup_target(test_par_conv)
target_link_libraries(test_par_conv fmt::fmt prm_handler parabolic_lib)

##############################################
#              My Executables                #
##############################################

##############################################
#                    Labs                    #
##############################################
#include(labs/L01/CMakeLists.txt)
#include(labs/L02/CMakeLists.txt)
#include(labs/L03/CMakeLists.txt)
#include(labs/L04/CMakeLists.txt)
#include(labs/L05/CMakeLists.txt)
#include(labs/L06/CMakeLists.txt)
#include(labs/L07/CMakeLists.txt)
#include(labs/L08/CMakeLists.txt)