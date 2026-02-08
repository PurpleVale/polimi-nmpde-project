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
add_executable(test1 test/test1.cpp)
deal_ii_setup_target(test1)
target_link_libraries(test1 fmt::fmt)

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