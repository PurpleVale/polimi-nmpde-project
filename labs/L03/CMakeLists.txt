add_executable(lab03ex01 labs/L03/exercise-01.cpp labs/L03/DiffusionReaction.cpp)
deal_ii_setup_target(lab03ex01)

add_executable(lab03ex01p labs/L03/exercise-01-parallel.cpp labs/L03/DiffusionReaction-parallel.cpp)
deal_ii_setup_target(lab03ex01p)
