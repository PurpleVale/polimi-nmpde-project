#include <Generics.hpp>
#include <PDEParamsHandler.hpp>
#include <ParabolicSolver.hpp>

int main(int argc , char** argv ) {

    dealii::Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

    unsigned int mpi_rank(dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD));
    dealii::ConditionalOStream pcout(std::cout, mpi_rank == 0);

    ParabolicPDE::ParabolicSolver<3> solver;
    solver.print_parameters("./text/parameters/parabolic_solver.prm");
    solver.print_editable_parameters("./text/parameters/parabolic_solver.xml");
    pcout << "================ " << fmt::format("{:^45}","STARTING TEST") << " ================" <<std::endl;
    double total_time = solver.run("./text/parameters/test_params/parabolic_solver_test.xml");
    pcout << fmt::format("{:<79}",fmt::format("Total time for solve: {:0.3f}s",total_time)) <<std::endl;
    pcout << "================ " << fmt::format("{:^45}","TEST DONE") << " ================" <<std::endl;
}