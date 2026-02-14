#ifndef NNPDE_STUDY_PDEPARAMSHANDLER_HPP
#define NNPDE_STUDY_PDEPARAMSHANDLER_HPP

#include <Generics.hpp>

namespace EllipticPDE {
    using namespace dealii;

    template <int dim>
    class EllipticParamHandler {

    public:
        using String = std::string;
        using ConstantMap = std::map<String, double>;
        using BoundaryIds  = types::boundary_id;

        EllipticParamHandler() :
        diffusion_c(1),
        advection_c(dim),
        reaction_c(1),
        force_term(1),
        mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)),
        mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)),
        pcout(std::cout, mpi_rank == 0) {

        }
        ~EllipticParamHandler() = default;

        void declare_parameters();
        void init(const String &filename);
        void print_parameters(const String &filename);
        void print_editable_parameters(const String &filename);

        String mesh_filename,output_filename;
        ParameterHandler prm;
        FunctionParser<dim> diffusion_c,
                            advection_c,
                            reaction_c,
                            force_term;
        unsigned int max_iters,poly_deg;
        String fe_type, preconditioner;
        std::vector<BoundaryIds> dirichlet_bc_tags, neumann_bc_tags;
        std::vector<std::unique_ptr<FunctionParser<dim>>> dirichlet_bc, neumann_bc;
        double epsilon;
        double residual_tolerance;
        bool symmetric_solver;
        bool param_initialized = false;
        bool initialized = false;

        private:
        void print_parameters_as(const String &filename,ParameterHandler::OutputStyle style);

        String variables = FunctionParser<dim>::default_variable_names();
        ConstantMap constants;

        const unsigned int mpi_size;
        const unsigned int mpi_rank;
        ConditionalOStream pcout;

    };

    template <int dim>
    class DNParamHandler {
    public:
        using String = std::string;
        using ConstantMap = std::map<String, double>;
        using BoundaryIds  = types::boundary_id;

        DNParamHandler() :
        mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)),
        mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)),
        pcout(std::cout, mpi_rank == 0)
        {};
        ~DNParamHandler() = default;

        void declare_parameters();
        void init(const String &filename);
        void print_parameters(const String &filename);
        void print_editable_parameters(const String &filename);

        String dirichlet_prm_file,neumann_prm_file;
        double relaxation;
        unsigned int max_iters;
        double tolerance;

        ParameterHandler prm;
        bool param_initialized = false;
        bool initialized = false;

        private:
        void print_parameters_as(const String &filename,ParameterHandler::OutputStyle style);

        const unsigned int mpi_size;
        const unsigned int mpi_rank;
        ConditionalOStream pcout;
    };

}
template class EllipticPDE::EllipticParamHandler<1>;
template class EllipticPDE::EllipticParamHandler<2>;
template class EllipticPDE::EllipticParamHandler<3>;
template class EllipticPDE::DNParamHandler<1>;
template class EllipticPDE::DNParamHandler<2>;
template class EllipticPDE::DNParamHandler<3>;


namespace ParabolicPDE {
    using namespace dealii;

    template <int dim>
    class ParabolicParamHandler {

    public:
        using String = std::string;
        using ConstantMap = std::map<String, double>;
        using BoundaryIds  = types::boundary_id;

        ParabolicParamHandler() :
        diffusion_c(1),
        advection_c(dim),
        reaction_c(1),
        force_term(1,0.0),
        initial_state(1),
        mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)),
        mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)),
        pcout(std::cout, mpi_rank == 0) {

        }

        ~ParabolicParamHandler() = default;

        void declare_parameters();
        void init(const String &filename);
        void print_parameters(const String &filename);
        void print_editable_parameters(const String &filename);

        String mesh_filename,output_filename;
        ParameterHandler prm;
        FunctionParser<dim> diffusion_c,
                            advection_c,
                            reaction_c,
                            force_term,
                            initial_state;
        unsigned int max_iters,poly_deg;
        String fe_type, preconditioner;
        std::vector<BoundaryIds> dirichlet_bc_tags, neumann_bc_tags;
        std::vector<std::unique_ptr<FunctionParser<dim>>> dirichlet_bc, neumann_bc;
        double epsilon;
        double residual_tolerance;
        double theta,time_end,time_step;
        bool symmetric_solver;
        bool param_initialized = false;
        bool initialized = false;

        private:
        void print_parameters_as(const String &filename,ParameterHandler::OutputStyle style);

        String variables = FunctionParser<dim>::default_variable_names() + ",t";
        String variables_no_time = FunctionParser<dim>::default_variable_names();
        ConstantMap constants;

        const unsigned int mpi_size;
        const unsigned int mpi_rank;
        ConditionalOStream pcout;

    };

}
template class ParabolicPDE::ParabolicParamHandler<1>;
template class ParabolicPDE::ParabolicParamHandler<2>;
template class ParabolicPDE::ParabolicParamHandler<3>;


namespace NavierStokesPDE {
    using namespace dealii;

    template <int dim>
    class StokesParamHandler {

    public:
        using String = std::string;
        using ConstantMap = std::map<String, double>;
        using BoundaryIds  = types::boundary_id;

        StokesParamHandler() :
        mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)),
        mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)),
        pcout(std::cout, mpi_rank == 0) {

        }
        ~StokesParamHandler() = default;

        void declare_parameters();
        void init(const String &filename);
        void print_parameters(const String &filename);
        void print_editable_parameters(const String &filename);

        String mesh_filename,output_filename;
        ParameterHandler prm;
        unsigned int max_iters,poly_deg_vel,poly_deg_pressure;
        double nu;  // [m^2/s]
        double sigma;
        double epsilon;
        bool param_initialized = false;
        bool initialized = false;

        private:
        void print_parameters_as(const String &filename,ParameterHandler::OutputStyle style);

        String variables = FunctionParser<dim>::default_variable_names();
        ConstantMap constants;

        const unsigned int mpi_size;
        const unsigned int mpi_rank;
        ConditionalOStream pcout;

    };


    template <int dim>
    class UnsteadyStokesParamHandler {

    public:
        using String = std::string;
        using ConstantMap = std::map<String, double>;
        using BoundaryIds  = types::boundary_id;

        UnsteadyStokesParamHandler() :
        mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)),
        mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)),
        pcout(std::cout, mpi_rank == 0) {

        }
        ~UnsteadyStokesParamHandler() = default;

        void declare_parameters();
        void init(const String &filename);
        void print_parameters(const String &filename);
        void print_editable_parameters(const String &filename);

        String mesh_filename,output_filename;
        ParameterHandler prm;
        unsigned int max_iters,poly_deg_vel,poly_deg_pressure;
        double nu;  // [m^2/s]
        double sigma;
        double epsilon;
        double T,dt,theta;
        bool param_initialized = false;
        bool initialized = false;

        private:
        void print_parameters_as(const String &filename,ParameterHandler::OutputStyle style);

        String variables = FunctionParser<dim>::default_variable_names();
        ConstantMap constants;

        const unsigned int mpi_size;
        const unsigned int mpi_rank;
        ConditionalOStream pcout;

    };

}
template class NavierStokesPDE::StokesParamHandler<1>;
template class NavierStokesPDE::StokesParamHandler<2>;
template class NavierStokesPDE::StokesParamHandler<3>;
template class NavierStokesPDE::UnsteadyStokesParamHandler<1>;
template class NavierStokesPDE::UnsteadyStokesParamHandler<2>;
template class NavierStokesPDE::UnsteadyStokesParamHandler<3>;

#endif //NNPDE_STUDY_PDEPARAMSHANDLER_HPP