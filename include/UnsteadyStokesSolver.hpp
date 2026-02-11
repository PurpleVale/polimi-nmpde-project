#ifndef NNPDE_STUDY_ELLIPTIC_SOLVER_HPP
#define NNPDE_STUDY_ELLIPTIC_SOLVER_HPP

#include <Generics.hpp>
#include <PDEParamsHandler.hpp>
#include <NavierStokesHelpers.hpp>

namespace NavierStokesPDE{
    using namespace dealii;

    template<int dim>
    class UnsteadyStokesSolver : public UnsteadyStokesParamHandler<dim>{

    public:
        using String = std::string;
        using BoundaryMap =  std::map<types::global_dof_index, double>;
        using BoundaryIds  = types::boundary_id;
        using BoundaryFunctionMap = std::map<types::boundary_id, const Function<dim> *>;
        using DTriangulation = parallel::fullydistributed::Triangulation<dim>;
        using BlkVector = TrilinosWrappers::MPI::BlockVector;
        using BlkSpMtx = TrilinosWrappers::BlockSparseMatrix;
        using DVector = TrilinosWrappers::MPI::Vector;
        using InterpretationVector = std::vector<DataComponentInterpretation::DataComponentInterpretation>;

        UnsteadyStokesSolver() :
            UnsteadyStokesParamHandler<dim>(),
            mesh(MPI_COMM_WORLD),
            force(dim),
            initial_sol(dim),
            d_bc_0(),
            d_bc_1(),
            n_bc_0(1),
            n_bc_1(1),
            mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)),
            mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)),
            pcout(std::cout, mpi_rank == 0)
         {}

        double setup(String parameter_filename);
        double assemble();
        double solve();
        void output() const;
        double run(const String &parameter_filename);
        // double compare_solution(
        //     const VectorTools::NormType &norm_type,
        //     const Function<dim> &exact
        // ) const;

    protected:

        DTriangulation mesh;

        std::unique_ptr<FiniteElement<dim>> fe;
        std::unique_ptr<Quadrature<dim>> quadrature;
        std::unique_ptr<Quadrature<dim-1>> boundary_quadrature;
        DoFHandler<dim> dof_handler;

        BlkSpMtx stiff;
        BlkSpMtx mass_pressure;
        BlkVector rhs;
        BlkVector sol;
        BlkVector owned_sol;

        IndexSet locally_owned_dofs;
        IndexSet locally_relevant_dofs;
        std::vector<IndexSet> block_owned_dofs;
        std::vector<IndexSet> block_relevant_dofs;

        StokesFunctions::Force<dim> force;
        StokesFunctions::InitialSolution<dim> initial_sol;
        StokesFunctions::DirichletCondition0<dim> d_bc_0;
        StokesFunctions::DirichletCondition1<dim> d_bc_1;
        StokesFunctions::NeumannCondition0<dim> n_bc_0;
        StokesFunctions::NeumannCondition1<dim> n_bc_1;

        double time;
        unsigned int time_step;

        const unsigned int mpi_size;
        const unsigned int mpi_rank;
        ConditionalOStream pcout;
    };

}
template class NavierStokesPDE::UnsteadyStokesSolver<1>;
template class NavierStokesPDE::UnsteadyStokesSolver<2>;
template class NavierStokesPDE::UnsteadyStokesSolver<3>;

#endif //NNPDE_STUDY_ELLIPTIC_SOLVER_HPP
