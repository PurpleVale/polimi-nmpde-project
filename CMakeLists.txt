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

##############################################
#                 My Tests                   #
##############################################
add_executable(test_ell test/test_elliptic.cpp)
deal_ii_setup_target(test_ell)
target_link_libraries(test_ell fmt::fmt)

add_executable(test_par test/test_parabolic.cpp)
deal_ii_setup_target(test_par)
target_link_libraries(test_par fmt::fmt)

add_executable(test_s_ell test/test_simple_elliptic.cpp)
deal_ii_setup_target(test_s_ell)
target_link_libraries(test_s_ell fmt::fmt)

add_executable(test_s_par test/test_simple_parabolic.cpp)
deal_ii_setup_target(test_s_par)
target_link_libraries(test_s_par fmt::fmt)
##############################################
#              My Executables                #
##############################################

##############################################
#                    Labs                    #
##############################################
include(labs/L01/CMakeLists.txt)
include(labs/L02/CMakeLists.txt)
include(labs/L03/CMakeLists.txt)
include(labs/L04/CMakeLists.txt)
include(labs/L05/CMakeLists.txt)
include(labs/L06/CMakeLists.txt)
include(labs/L07/CMakeLists.txt)
include(labs/L08/CMakeLists.txt)