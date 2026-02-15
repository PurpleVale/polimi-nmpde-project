#include <Generics.hpp>
#include <PDEParamsHandler.hpp>
#include <SimpleGridParabolicSolver.hpp>

class Exact : public dealii::Function<1> {

    public:
    Exact() {};

     virtual double value(const dealii::Point<1> &p,const unsigned int = 0) const override {
         return std::sin(M_PI*p[0]/2);
     }
    virtual dealii::Tensor<1, 1> gradient(const dealii::Point<1> &p, const unsigned int = 0) const override {
         dealii::Tensor<1, 1> grad;
         grad[0] = (M_PI/2) * std::cos(M_PI*p[0]/2);
         return grad;
     }
};

int main(int argc , char** argv ) {

    std::vector<double> dts = {0.1,0.05,0.025,0.0125};
    dealii::Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);
    std::ofstream outfile("./text/outputs/convergence.csv");
    outfile << "dt,eL2,eH1" << std::endl;

    unsigned int mpi_rank(dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD));
    dealii::ConditionalOStream pcout(std::cout, mpi_rank == 0);

    Exact exact;
    ParabolicPDE::SimpleGridParabolicSolver solver(10);
    dealii::ConvergenceTable table;

    solver.init("./text/parameters/test_params/parabolic_solver_test.xml");
    for (auto dt : dts) {
        solver.time_step = dt;
        solver.run("./text/parameters/test_params/parabolic_solver_test.xml");
        std::cout << "dt=" << solver.time_step << std::endl;
        auto L2 = solver.compare_solution(dealii::VectorTools::L2_norm,exact);
        auto H1 = solver.compare_solution(dealii::VectorTools::H1_norm,exact);

        table.add_value("dt",dt);
        table.add_value("L2",L2);
        table.add_value("H1",H1);

        outfile << dt << "," << L2 << "," << H1 << std::endl;
    }

    table.evaluate_all_convergence_rates(dealii::ConvergenceTable::reduction_rate_log2);
    table.set_scientific("L2", true);
    table.set_scientific("H1", true);

    table.write_text(std::cout);
}
