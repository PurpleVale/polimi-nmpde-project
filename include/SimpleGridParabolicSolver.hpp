#ifndef NNPDE_STUDY_ELLIPTIC_SOLVER_HPP
#define NNPDE_STUDY_ELLIPTIC_SOLVER_HPP

#include <Generics.hpp>
#include <PDEParamsHandler.hpp>

namespace ParabolicPDE{
    using namespace dealii;

    constexpr auto dim = 1;
    class SimpleGridParabolicSolver : public ParabolicParamHandler<dim>{

    public:
        using String = std::string;
        using BoundaryMap =  std::map<types::global_dof_index, double>;
        using BoundaryIds  = types::boundary_id;
        using BoundaryFunctionMap = std::map<types::boundary_id, const Function<dim> *>;

        SimpleGridParabolicSolver(unsigned int N_elm = 10) :
            ParabolicParamHandler<dim>(),
            mesh(MPI_COMM_WORLD),
            N_elm(N_elm),
            mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)),
            mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)),
            pcout(std::cout, mpi_rank == 0)
         {}

        double setup(const String parameter_filename);
        double assemble();
        double assemble_matrix();
        double assemble_rhs();
        double solve();
        void output() const;
        double run(const String &parameter_filename);
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
        TrilinosWrappers::MPI::Vector owned_sol;

        IndexSet locally_owned_dofs;
        IndexSet locally_relevant_dofs;

        double curr_time;
        unsigned int curr_time_step;
        unsigned int N_elm;

        const unsigned int mpi_size;
        const unsigned int mpi_rank;
        ConditionalOStream pcout;

    };

}

#endif //NNPDE_STUDY_ELLIPTIC_SOLVER_HPP
