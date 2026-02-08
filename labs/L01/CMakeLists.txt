add_library(Poisson1D STATIC labs/L01/Poisson1D.cpp)
deal_ii_setup_target(Poisson1D)

add_executable(lab01ex01 labs/L01/exercise-01.cpp)
target_link_libraries(lab01ex01 Poisson1D)
deal_ii_setup_target(lab01ex01)

add_executable(lab01ex02 labs/L01/exercise-02.cpp)
target_link_libraries(lab01ex02 Poisson1D)
deal_ii_setup_target(lab01ex02)
