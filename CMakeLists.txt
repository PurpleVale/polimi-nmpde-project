##############################################
#              Project Setup                 #
##############################################
cmake_minimum_required(VERSION 3.24.0)
project(NnPDE_Project LANGUAGES CXX C)

include(./cmake-common.cmake)
include_directories(include)

##############################################
#              My Libraries                  #
##############################################


##############################################
#                 My Tests                   #
##############################################
add_executable(test1 test/test1.cpp)
deal_ii_setup_target(test1)

##############################################
#              My Executables                #
##############################################
#add_executable()
#deal_ii_setup_target()
