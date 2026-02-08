#ifndef NMPDE_PROJECT_GENERICS_HPP
#define NMPDE_PROJECT_GENERICS_HPP

//**********************************************/
/***             General Includes            ***/
/***********************************************/
//this is bad programming but it's fine
#include <deal.II/base/timer.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/function_parser.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <fmt/format.h>

#if !defined(BUILD_TYPE_DEBUG)
    #define LOG_VAR(name,val) ;
    #define LOG(val) ;
    #define LOG_TITLE(title) ;
    #define LOG_ANY(formt,arg) ;
#else
    #define GREEN   "\033[32m"
    #define YELLOW  "\033[33m"
    #define BLUE  "\033[34m"
    #define CYAN    "\033[36m"
    #define RESET   "\033[0m"

    #define FILE_NAME fmt::format("{:<21}",std::filesystem::path(__FILE__).filename().string())
    #define LOG_PREFIX  CYAN "LOG " << BLUE << FILE_NAME << CYAN << " at line " << BLUE << fmt::format("{:>04d}",__LINE__) << CYAN << ":\t [ "
    #define LOG_SUFFIX  " ]" << RESET
    #define LOG_TITLE(title) pcout << LOG_PREFIX << CYAN \
        << "================ " << fmt::format("{:^45}",title) << " ================" \
        << LOG_SUFFIX <<std::endl;
    #define LOG_VAR(name,val) pcout << LOG_PREFIX << GREEN << fmt::format("{:<89}",fmt::format("{}{} == {}{}",name, CYAN, YELLOW, val)) << CYAN << LOG_SUFFIX <<std::endl;
    #define LOG(val) pcout << LOG_PREFIX << YELLOW << fmt::format("{:<89}",val) << CYAN << LOG_SUFFIX <<std::endl;
    #define LOG_ANY(formt,arg) pcout << LOG_PREFIX << YELLOW << fmt::format("{:<79}",fmt::format(formt,arg)) << CYAN << LOG_SUFFIX <<std::endl;
#endif

#endif //NMPDE_PROJECT_GENERICS_HPP