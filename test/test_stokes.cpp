#include <Generics.hpp>
#include <PDEParamsHandler.hpp>
#include <StokesSolver.hpp>

int main(int argc , char** argv ) {

    dealii::Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

    unsigned int mpi_rank(dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD));
    dealii::ConditionalOStream pcout(std::cout, mpi_rank == 0);

    NavierStokesPDE::StokesSolver<3> solver;
    solver.print_parameters("./text/parameters/stokes_solver.prm");
    solver.print_editable_parameters("./text/parameters/stokes_solver.xml");
    pcout << "================ " << fmt::format("{:^45}","STARTING TEST") << " ================" <<std::endl;
    double setup_time = solver.setup("./text/parameters/test_params/stokes_solver_test.xml");
    double assemble_time = solver.assemble();
    double solve_time = solver.solve();
    solver.output();
    pcout << fmt::format("{:<79}",fmt::format("Setup time: {:0.3f}s",setup_time)) <<std::endl;
    pcout << fmt::format("{:<79}",fmt::format("Assemble time: {:0.3f}s",assemble_time)) <<std::endl;
    pcout << fmt::format("{:<79}",fmt::format("Solve time: {:0.3f}s",solve_time)) <<std::endl;
    pcout << "================ " << fmt::format("{:^45}","TEST DONE") << " ================" <<std::endl;
}