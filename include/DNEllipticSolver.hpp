#ifndef NNPDE_STUDY_ELLIPTIC_SOLVER_HPP
#define NNPDE_STUDY_ELLIPTIC_SOLVER_HPP

#include <Generics.hpp>
#include <PDEParamsHandler.hpp>

namespace EllipticPDE{
    using namespace dealii;

    template<int dim>
    class DNEllipticSolver : public EllipticParamHandler<dim>{

    public:
        using String = std::string;
        using IndexMap = std::map<types::global_dof_index, types::global_dof_index>;
        using PointMap = std::map<types::global_dof_index, Point<dim>>;
        using BoundaryMap =  std::map<types::global_dof_index, double>;
        using BoundaryIds  = types::boundary_id;
        using BoundaryFunctionMap = std::map<types::boundary_id, const Function<dim> *>;
        using DVector = TrilinosWrappers::MPI::Vector;
        using SpMtx = TrilinosWrappers::SparseMatrix;

        DNEllipticSolver(const unsigned int &subdomain_id) :
            EllipticParamHandler<dim>(),
            domain_id(subdomain_id),
            mesh(MPI_COMM_WORLD),
            mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)),
            mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)),
            pcout(std::cout, mpi_rank == 0)
         {}

        double setup(String parameter_filename);
        double assemble();
        double solve();
        void output(const unsigned int &iter) const;
        void apply_interface_dirichlet(const DNEllipticSolver &other);
        void apply_interface_neumann(DNEllipticSolver &other);
        const DVector& get_solution() const;
        void apply_relaxation(const DVector& old_solution, const double &lambda);
        IndexMap compute_interface_map(const DNEllipticSolver &other) const;
        double compare_solution(
            const VectorTools::NormType &norm_type,
            const Function<dim> &exact
        ) const;

    protected:

        unsigned int domain_id;

        parallel::fullydistributed::Triangulation<dim> mesh;

        PointMap support_points;

        std::unique_ptr<FiniteElement<dim>> fe;
        std::unique_ptr<Quadrature<dim>> quadrature;
        std::unique_ptr<Quadrature<dim-1>> boundary_quadrature;
        DoFHandler<dim> dof_handler;

        SpMtx stiff;
        DVector rhs;
        DVector sol;

        IndexSet locally_owned_dofs;

        const unsigned int mpi_size;
        const unsigned int mpi_rank;
        ConditionalOStream pcout;

    };

}
template class EllipticPDE::DNEllipticSolver<1>;
template class EllipticPDE::DNEllipticSolver<2>;
template class EllipticPDE::DNEllipticSolver<3>;

#endif //NNPDE_STUDY_ELLIPTIC_SOLVER_HPP
