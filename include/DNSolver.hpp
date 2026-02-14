#ifndef NNPDE_PROJECT_DNSOLVER_HPP
#define NNPDE_PROJECT_DNSOLVER_HPP

#include <Generics.hpp>
#include <DNEllipticSolver.hpp>

namespace EllipticPDE {

    template<int dim>
    class DNSolver : public DNParamHandler<dim> {
        public:
        using String = std::string;

        DNSolver() :
            DNParamHandler<dim>(),
            dirichlet_domain(0),
            neumann_domain(1),
            iter(0),
            mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)),
            mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)),
            pcout(std::cout, mpi_rank == 0)
        {};
        ~DNSolver() = default;

        DNEllipticSolver<dim> dirichlet_domain;
        DNEllipticSolver<dim> neumann_domain;

        unsigned int iter;

        double setup(String parameter_filename);
        double solve();

        protected:
        const unsigned int mpi_size;
        const unsigned int mpi_rank;
        ConditionalOStream pcout;

    };

}
template class EllipticPDE::DNSolver<1>;
template class EllipticPDE::DNSolver<2>;
template class EllipticPDE::DNSolver<3>;

#endif //NNPDE_PROJECT_DNSOLVER_HPP