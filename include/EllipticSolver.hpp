#ifndef NNPDE_STUDY_ELLIPTIC_SOLVER_HPP
#define NNPDE_STUDY_ELLIPTIC_SOLVER_HPP

#include <Generics.hpp>
#include <PDEParamsHandler.hpp>

namespace EllipticPDE{
    using namespace dealii;

    template<int dim>
    class EllipticSolver : public EllipticParamHandler<dim>{

    public:
        using String = std::string;
        using BoundaryMap =  std::map<types::global_dof_index, double>;
        using BoundaryIds  = types::boundary_id;
        using BoundaryFunctionMap = std::map<types::boundary_id, const Function<dim> *>;

        EllipticSolver() :
            EllipticParamHandler<dim>(),
            mesh(MPI_COMM_WORLD),
            mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)),
            mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)),
            pcout(std::cout, mpi_rank == 0)
         {}

        double setup(String parameter_filename);
        double assemble();
        double solve();
        void output() const;
        double compare_solution(
            const VectorTools::NormType &norm_type,
            const Function<dim> &exact
        ) const;

    protected:
        parallel::fullydistributed::Triangulation<dim> mesh;

        std::unique_ptr<FiniteElement<dim>> fe;
        std::unique_ptr<Quadrature<dim>> quadrature;
        std::unique_ptr<Quadrature<dim-1>> boundary_quadrature;
        DoFHandler<dim> dof_handler;

        TrilinosWrappers::SparseMatrix stiff;
        TrilinosWrappers::MPI::Vector rhs;
        TrilinosWrappers::MPI::Vector sol;

        IndexSet locally_owned_dofs;

        const unsigned int mpi_size;
        const unsigned int mpi_rank;
        ConditionalOStream pcout;

    };

}

template class EllipticPDE::EllipticSolver<1>;
template class EllipticPDE::EllipticSolver<2>;
template class EllipticPDE::EllipticSolver<3>;

#endif //NNPDE_STUDY_ELLIPTIC_SOLVER_HPP
