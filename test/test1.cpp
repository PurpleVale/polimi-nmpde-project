#include <Generics.hpp>
#include <PDEParamsHandler.hpp>
#include <EllipticSolver.hpp>

int main(int argc , char** argv ) {

    dealii::Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

    dealii::ConditionalOStream pcout(std::cout,MPI_COMM_WORLD);

    EllipticPDE::EllipticSolver<3> solver;
    solver.print_parameters("./text/parameters/elliptic_solver.prm");
    solver.print_editable_parameters("./text/parameters/elliptic_solver.xml");
    LOG_TITLE("STARTING TEST 1")
    double setup_time = solver.setup("./text/parameters/elliptic_solver_test1.xml");
    double assemble_time = solver.assemble();
    double solve_time = solver.solve();
    solver.output();
    LOG_ANY("Setup time: {:0.3f}s",setup_time)
    LOG_ANY("Assemble time: {:0.3f}s",assemble_time)
    LOG_ANY("Solve time: {:0.3f}s",solve_time)
    LOG_TITLE("TEST 1 DONE")
}