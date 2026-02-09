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
        dirichlet_bc(1),
        neumann_bc(1),
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
                            force_term,
                            dirichlet_bc,
                            neumann_bc;
        unsigned int max_iters,poly_deg;
        String quad_type, fe_type, preconditioner;
        std::vector<BoundaryIds> boundary_tags;
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

    template<int dim>
    void EllipticParamHandler<dim>::declare_parameters() {
        LOG_TITLE("Constructing Parameters Structure")
        prm.enter_subsection("Solver");
        {
            prm.declare_entry(
                "Max iters",
                "1000",
                Patterns::Integer(1,100000),
                "Max number of iters for linear solver"
            );
            prm.declare_entry(
                "Tolerance",
                "1e-10",
                Patterns::Double(),
                "Tolerance for the linear solver"
            );
            prm.declare_entry(
                "Residual tolerance",
                "1e-10",
                Patterns::Double(),
                "Tolerance for the residual of the linear solver"
            );
            prm.declare_entry(
                "Solver type",
                "CG",
                Patterns::Selection("CG|GMRES"),
                "The solver to use"
            );
            prm.declare_entry(
                "Preconditioner",
                "Identity",
                Patterns::Selection("Jacobi|Identity|ILU|ILUT|SSOR"),
                "The type of preconditioner default is none"
            );
        }
        prm.leave_subsection();
        prm.enter_subsection("Files");
        {
            prm.declare_entry(
                "Mesh file",
                "./mesh.vtk",
                Patterns::FileName(Patterns::FileName::input),
                "A .vtk file containing a mesh",
                true
            );
            prm.declare_entry(
                "Output file",
                "./output.vtk",
                Patterns::FileName(Patterns::FileName::output),
                "A .vtk file to save the solution to"
            );
        }
        prm.leave_subsection();
        prm.enter_subsection("Functions");
        {
            prm.declare_entry(
                "Diffusion",
                "1.0",
                Patterns::Anything(),
                "Diffusion coefficient"
            );
            prm.declare_entry(
                "Advection x",
                "0.0",
                Patterns::Anything(),
                "Advection x coefficient"
            );
            prm.declare_entry(
                "Advection y",
                "0.0",
                Patterns::Anything(),
                "Advection y coefficient"
            );
            prm.declare_entry(
                "Advection z",
                "0.0",
                Patterns::Anything(),
                "Advection z coefficient"
            );
            prm.declare_entry(
                "Reaction",
                "0.0",
                Patterns::Anything(),
                "Reaction coefficient"
            );
            prm.declare_entry(
                "Force",
                "0.0",
                Patterns::Anything(),
                "Force term"
            );
            prm.declare_entry(
                "Dirichlet BC",
                "0.0",
                Patterns::Anything(),
                "Dirichlet Boundary Condition"
            );
            prm.declare_entry(
                "Neumann BC",
                "0.0",
                Patterns::Anything(),
                "Neumann Boundary Condition"
            );
        }
        prm.leave_subsection();
        prm.enter_subsection("Finite Elements");
        {
            prm.declare_entry(
                "Polynomial degree",
                "2",
                Patterns::Integer(1,7),
                "The polynomial degree of the Elements"
            );
            prm.declare_entry(
                "Quadrature",
                "Gaussian",
                Patterns::Selection("Gaussian"),
                "The quadrature algorithm to be used"
            );
            prm.declare_entry(
                "Elements type",
                "Simplex",
                Patterns::Selection("Simplex"),
                "The finite elements to be used"
            );
            prm.declare_entry(
                "Dirichlet Tags",
                "0,1",
                Patterns::List(Patterns::Integer(0,9),0),
                "The Tag of the boundaries on which to apply the Dirichlet condition"
            );
        }
        prm.leave_subsection();
        param_initialized = true;
    }

    template<int dim>
    void EllipticParamHandler<dim>::init(const String &filename) {
        if (initialized) return;

        if (!param_initialized) {
            declare_parameters();
        }

        LOG_TITLE("Initializing Elliptic parameters");
        constants["pi"] = numbers::PI;

        String diffusion_s;
        std::vector<String> advection_s(dim);
        String reaction_s;
        String force_term_s;
        String dirichlet_bc_s;
        String neumann_bc_s;

        String tags;

        LOG_TITLE("Reading Parameter File")
        LOG_VAR("File Name", filename)
        prm.parse_input(filename);

        prm.enter_subsection("Solver");
        {
            LOG_TITLE("Getting Solver Parameters")
            max_iters = prm.get_integer("Max iters");
            epsilon = prm.get_double("Tolerance");
            residual_tolerance = prm.get_double("Residual tolerance");
            symmetric_solver = prm.get("Solver type") == "CG";
            preconditioner = prm.get("Preconditioner");
        }
        prm.leave_subsection();
        prm.enter_subsection("Files");
        {
            LOG_TITLE("Getting File Names")
            mesh_filename = prm.get("Mesh file");
            output_filename = prm.get("Output file");
        }
        prm.leave_subsection();
        prm.enter_subsection("Functions");
        {
            LOG_TITLE("Getting Functions")
            diffusion_s = prm.get("Diffusion");
            for (int d = 0; d < dim; d++) {
                String advection_tag;
                switch (d) {
                    case 0: advection_tag = "Advection x"; break;
                    case 1: advection_tag = "Advection y"; break;
                    case 2: advection_tag = "Advection z"; break;
                    default: ;
                }
                advection_s[d] = prm.get(advection_tag);
            }
            reaction_s = prm.get("Reaction");
            force_term_s = prm.get("Force");
            dirichlet_bc_s = prm.get("Dirichlet BC");
            neumann_bc_s = prm.get("Neumann BC");

            LOG_TITLE("Initializing Functions from strings")
            diffusion_c.initialize(variables,diffusion_s,constants);
            advection_c.initialize(variables,advection_s,constants);
            reaction_c.initialize(variables,reaction_s,constants);
            force_term.initialize(variables,force_term_s,constants);
            dirichlet_bc.initialize(variables,dirichlet_bc_s,constants);
            neumann_bc.initialize(variables,neumann_bc_s,constants);
        }
        prm.leave_subsection();
        prm.enter_subsection("Finite Elements");
        {
            LOG_TITLE("Getting Finite Elements Parameters")
            poly_deg = prm.get_integer("Polynomial degree");
            quad_type = prm.get("Quadrature");
            fe_type = prm.get("Elements type");
            tags = prm.get("Dirichlet Tags");
            for (auto tag : tags) {
                if (tag == ',') continue;
                boundary_tags.emplace_back(static_cast<BoundaryIds>(tag-'0'));
            }
        }
        prm.leave_subsection();

        LOG_TITLE("Parameters Read successful")
        LOG_VAR("Max iters",max_iters)
        LOG_VAR("Tolerance",epsilon)
        LOG_VAR("Residual tolerance",residual_tolerance)
        LOG_VAR("Solver", (symmetric_solver ? "CG":"GMRES") )
        LOG_VAR("Mesh file",mesh_filename)
        LOG_VAR("Output file",output_filename)
        LOG_VAR("Diffusion Coefficient",diffusion_s)
        for (int d = 0; d < dim; d++) {
            String advection_tag;
            switch (d) {
                case 0: advection_tag = "Advection x"; break;
                case 1: advection_tag = "Advection y"; break;
                case 2: advection_tag = "Advection z"; break;
                default: ;
            }
            LOG_VAR(advection_tag, advection_s[d]);
        }
        LOG_VAR("Reaction Coefficient",reaction_s)
        LOG_VAR("Force term", force_term_s)
        LOG_VAR("Dirichlet BC", dirichlet_bc_s)
        LOG_VAR("Neumann BC", neumann_bc_s)
        LOG_VAR("Polynomial Degree", poly_deg)
        LOG_VAR("Quadrature Type",quad_type)
        LOG_VAR("Finite Elements Type",fe_type)
        LOG_VAR("Tags for Dirichlet boundaries",tags)

        initialized = true;
    }

    template<int dim>
    void EllipticParamHandler<dim>::print_parameters(const String &filename) {
        print_parameters_as(filename,ParameterHandler::DefaultStyle);
    }

    template<int dim>
    void EllipticParamHandler<dim>::print_editable_parameters(const String &filename) {
        print_parameters_as(filename,ParameterHandler::XML);
    }

    template<int dim>
    void EllipticParamHandler<dim>::print_parameters_as(
        const String &filename,
        ParameterHandler::OutputStyle style
    ) {
        if (! param_initialized) {
            declare_parameters();
        }
        prm.print_parameters(filename,style);
    }

}


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
        dirichlet_bc(1,0.0),
        neumann_bc(1,0.0),
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
                            dirichlet_bc,
                            neumann_bc,
                            initial_state;
        unsigned int max_iters,poly_deg;
        String quad_type, fe_type, preconditioner;
        std::vector<BoundaryIds> boundary_tags;
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

    template<int dim>
    void ParabolicParamHandler<dim>::declare_parameters() {
        LOG_TITLE("Constructing Parameters Structure")
        prm.enter_subsection("Solver");
        {
            prm.declare_entry(
                "Max iters",
                "1000",
                Patterns::Integer(1,100000),
                "Max number of iters for linear solver"
            );
            prm.declare_entry(
                "Tolerance",
                "1e-10",
                Patterns::Double(),
                "Tolerance for the linear solver"
            );
            prm.declare_entry(
                "Residual tolerance",
                "1e-10",
                Patterns::Double(),
                "Tolerance for the residual of the linear solver"
            );
            prm.declare_entry(
                "Solver type",
                "CG",
                Patterns::Selection("CG|GMRES"),
                "The solver to use"
            );
            prm.declare_entry(
                "Preconditioner",
                "Identity",
                Patterns::Selection("Jacobi|Identity|ILU|ILUT|SSOR"),
                "The type of preconditioner default is none"
            );
        }
        prm.leave_subsection();
        prm.enter_subsection("Files");
        {
            prm.declare_entry(
                "Mesh file",
                "./mesh.vtk",
                Patterns::FileName(Patterns::FileName::input),
                "A .vtk file containing a mesh",
                true
            );
            prm.declare_entry(
                "Output file",
                "./output.vtk",
                Patterns::FileName(Patterns::FileName::output),
                "A .vtk file to save the solution to"
            );
        }
        prm.leave_subsection();
        prm.enter_subsection("Functions");
        {
            prm.declare_entry(
                "Diffusion",
                "1.0",
                Patterns::Anything(),
                "Diffusion coefficient"
            );
            prm.declare_entry(
                "Advection x",
                "0.0",
                Patterns::Anything(),
                "Advection x coefficient"
            );
            prm.declare_entry(
                "Advection y",
                "0.0",
                Patterns::Anything(),
                "Advection y coefficient"
            );
            prm.declare_entry(
                "Advection z",
                "0.0",
                Patterns::Anything(),
                "Advection z coefficient"
            );
            prm.declare_entry(
                "Reaction",
                "0.0",
                Patterns::Anything(),
                "Reaction coefficient"
            );
            prm.declare_entry(
                "Force",
                "0.0",
                Patterns::Anything(),
                "Force term"
            );
            prm.declare_entry(
                "Dirichlet BC",
                "0.0",
                Patterns::Anything(),
                "Dirichlet Boundary Condition"
            );
            prm.declare_entry(
                "Neumann BC",
                "0.0",
                Patterns::Anything(),
                "Neumann Boundary Condition"
            );
            prm.declare_entry(
                "Initial State",
                "0.0",
                Patterns::Anything(),
                "The initial state of the system"
            );
        }
        prm.leave_subsection();
        prm.enter_subsection("Finite Elements");
        {
            prm.declare_entry(
                "Polynomial degree",
                "2",
                Patterns::Integer(1,7),
                "The polynomial degree of the Elements"
            );
            prm.declare_entry(
                "Quadrature",
                "Gaussian",
                Patterns::Selection("Gaussian"),
                "The quadrature algorithm to be used"
            );
            prm.declare_entry(
                "Elements type",
                "Simplex",
                Patterns::Selection("Simplex"),
                "The finite elements to be used"
            );
            prm.declare_entry(
                "Dirichlet Tags",
                "0,1",
                Patterns::List(Patterns::Integer(0,9),0),
                "The Tag of the boundaries on which to apply the Dirichlet condition"
            );
        }
        prm.leave_subsection();
        prm.enter_subsection("Time Parameters");
        {
            prm.declare_entry(
                "Final time",
                "1.0",
                Patterns::Double(),
                "The T st the pde has domain Omega times (0,T)",
                true
            );
            prm.declare_entry(
                "Time step",
                "0.05",
                Patterns::Double(),
                "The time step dt to use for the theta method"
            );
            prm.declare_entry(
                "Theta",
                "0.5",
                Patterns::Double(0.0,1.0),
                "The theta parameter for the theta method"
            );
        }
        prm.leave_subsection();
        param_initialized = true;
    }

    template<int dim>
    void ParabolicParamHandler<dim>::init(const String &filename) {
        if (initialized) return;

        if (!param_initialized) {
            declare_parameters();
        }

        LOG_TITLE("Initializing Parabolic parameters");
        constants["pi"] = numbers::PI;

        String diffusion_s;
        std::vector<String> advection_s(dim);
        String reaction_s;
        String force_term_s;
        String dirichlet_bc_s;
        String neumann_bc_s;
        String init_state_s;

        String tags;

        LOG_TITLE("Reading Parameter File")
        LOG_VAR("File Name", filename)
        prm.parse_input(filename);

        prm.enter_subsection("Solver");
        {
            LOG_TITLE("Getting Solver Parameters")
            max_iters = prm.get_integer("Max iters");
            epsilon = prm.get_double("Tolerance");
            residual_tolerance = prm.get_double("Residual tolerance");
            symmetric_solver = prm.get("Solver type") == "CG";
            preconditioner = prm.get("Preconditioner");
        }
        prm.leave_subsection();
        prm.enter_subsection("Files");
        {
            LOG_TITLE("Getting File Names")
            mesh_filename = prm.get("Mesh file");
            output_filename = prm.get("Output file");
        }
        prm.leave_subsection();
        prm.enter_subsection("Functions");
        {
            LOG_TITLE("Getting Functions")
            diffusion_s = prm.get("Diffusion");
            for (int d = 0; d < dim; d++) {
                String advection_tag;
                switch (d) {
                    case 0: advection_tag = "Advection x"; break;
                    case 1: advection_tag = "Advection y"; break;
                    case 2: advection_tag = "Advection z"; break;
                    default: ;
                }
                advection_s[d] = prm.get(advection_tag);
            }
            reaction_s = prm.get("Reaction");
            force_term_s = prm.get("Force");
            dirichlet_bc_s = prm.get("Dirichlet BC");
            neumann_bc_s = prm.get("Neumann BC");
            init_state_s = prm.get("Initial State");

            LOG_TITLE("Initializing Functions from strings")
            diffusion_c.initialize(variables_no_time,diffusion_s,constants);
            advection_c.initialize(variables_no_time,advection_s,constants);
            reaction_c.initialize(variables_no_time,reaction_s,constants);
            force_term.initialize(variables,force_term_s,constants,true);
            dirichlet_bc.initialize(variables,dirichlet_bc_s,constants,true);
            neumann_bc.initialize(variables,neumann_bc_s,constants,true);
            initial_state.initialize(variables_no_time,init_state_s,constants);
        }
        prm.leave_subsection();
        prm.enter_subsection("Finite Elements");
        {
            LOG_TITLE("Getting Finite Elements Parameters")
            poly_deg = prm.get_integer("Polynomial degree");
            quad_type = prm.get("Quadrature");
            fe_type = prm.get("Elements type");
            tags = prm.get("Dirichlet Tags");
            for (auto tag : tags) {
                if (tag == ',') continue;
                boundary_tags.emplace_back(static_cast<BoundaryIds>(tag-'0'));
            }
        }
        prm.leave_subsection();
        prm.enter_subsection("Time Parameters");
        {
            LOG_TITLE("Getting Time Parameters")
            time_end = prm.get_double("Final time");
            time_step = prm.get_double("Time step");
            theta = prm.get_double("Theta");
        }
        prm.leave_subsection();

        LOG_TITLE("Parameters Read successful")
        LOG_VAR("Max iters",max_iters)
        LOG_VAR("Tolerance",epsilon)
        LOG_VAR("Residual tolerance",residual_tolerance)
        LOG_VAR("Solver", (symmetric_solver ? "CG":"GMRES") )
        LOG_VAR("Mesh file",mesh_filename)
        LOG_VAR("Output file",output_filename)
        LOG_VAR("Diffusion Coefficient",diffusion_s)
        for (int d = 0; d < dim; d++) {
            String advection_tag;
            switch (d) {
                case 0: advection_tag = "Advection x"; break;
                case 1: advection_tag = "Advection y"; break;
                case 2: advection_tag = "Advection z"; break;
                default: ;
            }
            LOG_VAR(advection_tag, advection_s[d]);
        }
        LOG_VAR("Reaction Coefficient",reaction_s)
        LOG_VAR("Force term", force_term_s)
        LOG_VAR("Dirichlet BC", dirichlet_bc_s)
        LOG_VAR("Neumann BC", neumann_bc_s)
        LOG_VAR("Polynomial Degree", poly_deg)
        LOG_VAR("Quadrature Type",quad_type)
        LOG_VAR("Finite Elements Type",fe_type)
        LOG_VAR("Tags for Dirichlet boundaries",tags)
        LOG_VAR("Time limit of the PDE",time_end)
        LOG_VAR("Time step",time_step)
        LOG_VAR("Theta Method Parameter",theta)

        initialized = true;
    }

    template<int dim>
    void ParabolicParamHandler<dim>::print_parameters(const String &filename) {
        print_parameters_as(filename,ParameterHandler::DefaultStyle);
    }

    template<int dim>
    void ParabolicParamHandler<dim>::print_editable_parameters(const String &filename) {
        print_parameters_as(filename,ParameterHandler::XML);
    }

    template<int dim>
    void ParabolicParamHandler<dim>::print_parameters_as(
        const String &filename,
        ParameterHandler::OutputStyle style
    ) {
        if (! param_initialized) {
            declare_parameters();
        }
        prm.print_parameters(filename,style);
    }

}



#endif //NNPDE_STUDY_PDEPARAMSHANDLER_HPP