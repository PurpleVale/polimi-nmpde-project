#include <Generics.hpp>
#include <PDEParamsHandler.hpp>
#include <ParabolicSolver.hpp>

class Exact : public dealii::Function<2> {

    public:
    Exact() {};

     virtual double value(const dealii::Point<2> &p,const unsigned int = 0) const override {
         return 0;
     }
    virtual dealii::Tensor<1, 2> gradient(const dealii::Point<2> &p, const unsigned int = 0) const override {
         dealii::Tensor<1, 2> grad;
         grad[0] = 0;
         grad[1] = 0;
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
    ParabolicPDE::ParabolicSolver<2> solver;
    dealii::ConvergenceTable table;

    solver.init("./text/parameters/parabolic.prm");
    for (auto dt : dts) {
        solver.time_step = dt;
        solver.run("./text/parameters/parabolic.prm");
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
