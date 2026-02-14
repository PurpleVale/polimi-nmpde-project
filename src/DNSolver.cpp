#include <Generics.hpp>
#include <PDEParamsHandler.hpp>
#include <DNSolver.hpp>

namespace EllipticPDE {
    template<int dim>
    double DNSolver<dim>::setup(String parameter_filename) {
        this->init(parameter_filename);

        double total_time = 0.0;
        LOG_TITLE("Initializing Dirichlet Domain")
        total_time += dirichlet_domain.setup(this->dirichlet_prm_file);
        LOG_TITLE("Initializing Neumann Domain")
        total_time += neumann_domain.setup(this->neumann_prm_file);
        return total_time;
    }

    template<int dim>
    double DNSolver<dim>::solve() {
        const bool print_iter = true;
        double total_time = 0.0;

        double increment_norm = this->tolerance + 1;
        iter = 0;

        while (iter < this->max_iters && increment_norm > this->tolerance) {

            if (print_iter)
                pcout << fmt::format("{:<90}",fmt::format("{} == {}","Starting iter", iter)) << std::endl;

            auto curr_sol_neumann = neumann_domain.get_solution();

            total_time += dirichlet_domain.assemble();
            dirichlet_domain.apply_interface_dirichlet(neumann_domain);
            total_time += dirichlet_domain.solve();

            total_time += neumann_domain.assemble();
            neumann_domain.apply_interface_neumann(dirichlet_domain);
            total_time += neumann_domain.solve();

            neumann_domain.apply_relaxation(curr_sol_neumann,this->relaxation);

            // curr_sol_neumann = curr_sol_neumann - neumann_domain.get_solution()
            curr_sol_neumann.add(-1.0,neumann_domain.get_solution());
            increment_norm = curr_sol_neumann.l2_norm();

            dirichlet_domain.output(iter);
            neumann_domain.output(iter);

            ++iter;
        }
        return total_time;
    }
}