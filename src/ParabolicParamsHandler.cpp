#include <Generics.hpp>
#include <PDEParamsHandler.hpp>

namespace ParabolicPDE {
    using namespace dealii;

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
                Patterns::List(Patterns::Anything(),0),
                "Dirichlet Boundary Condition"
            );
            prm.declare_entry(
                "Neumann BC",
                "0.0",
                Patterns::List(Patterns::Anything(),0),
                "Neumann Boundary Condition"
            );
            prm.declare_entry(
                "Dirichlet Tags",
                "0,1",
                Patterns::List(Patterns::Integer(0,9),0),
                "The Tag of the boundaries on which to apply the Dirichlet condition"
            );
            prm.declare_entry(
                "Neumann Tags",
                "2,3",
                Patterns::List(Patterns::Integer(0,9),0),
                "The Tag of the boundaries on which to apply the Neumann condition"
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
                "Elements type",
                "Simplex",
                Patterns::Selection("Simplex|Q"),
                "The finite elements to be used"
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
        String init_state_s;
        std::vector<String> dirichlet_bc_s;
        std::vector<String> neumann_bc_s;

        std::vector<BoundaryIds> tags_d;
        std::vector<BoundaryIds> tags_n;

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

            String dirichlet_bc_curr, neumann_bc_curr;

            std::stringstream dirichlet_bc_ss(prm.get("Dirichlet BC"));
            std::stringstream dirichlet_tag_ss(prm.get("Dirichlet Tags"));
            while (std::getline(dirichlet_bc_ss, dirichlet_bc_curr, ',')) {
                dirichlet_bc_s.emplace_back(dirichlet_bc_curr);
            }
            while (std::getline(dirichlet_tag_ss, dirichlet_bc_curr, ',')) {
                dirichlet_bc_tags.emplace_back(static_cast<BoundaryIds>(stoi(dirichlet_bc_curr)));
            }
            for (long unsigned int i = 0; i < dirichlet_bc_tags.size(); i++) {
                dirichlet_bc.push_back(std::make_unique<FunctionParser<dim>>(1,0.0));
                dirichlet_bc[i]->initialize(variables,dirichlet_bc_s[i],constants,true);
            }
            std::stringstream neumann_bc_ss(prm.get("Neumann BC"));
            std::stringstream neumann_tag_ss(prm.get("Neumann Tags"));
            while (std::getline(neumann_bc_ss, neumann_bc_curr, ',')) {
                neumann_bc_s.emplace_back(neumann_bc_curr);
            }
            while (std::getline(neumann_tag_ss, neumann_bc_curr, ',')) {
                neumann_bc_tags.emplace_back(static_cast<BoundaryIds>(stoi(neumann_bc_curr)));
            }
            for (long unsigned int i = 0; i < neumann_bc_tags.size(); i++) {
                neumann_bc.push_back(std::make_unique<FunctionParser<dim>>(1,0.0));
                neumann_bc[i]->initialize(variables,neumann_bc_s[i],constants,true);
            }

            init_state_s = prm.get("Initial State");

            LOG_TITLE("Initializing Functions from strings")
            diffusion_c.initialize(variables_no_time,diffusion_s,constants);
            advection_c.initialize(variables_no_time,advection_s,constants);
            reaction_c.initialize(variables_no_time,reaction_s,constants);
            force_term.initialize(variables,force_term_s,constants,true);
            initial_state.initialize(variables_no_time,init_state_s,constants);
        }
        prm.leave_subsection();
        prm.enter_subsection("Finite Elements");
        {
            LOG_TITLE("Getting Finite Elements Parameters")
            poly_deg = prm.get_integer("Polynomial degree");
            fe_type = prm.get("Elements type");
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
        for (long unsigned int i = 0; i < dirichlet_bc_tags.size(); i++) {
            LOG(fmt::format("Dirichlet Boundary {}, Function: {}",dirichlet_bc_tags[i],dirichlet_bc_s[i]));
        }
        for (long unsigned int i = 0; i < neumann_bc_tags.size(); i++) {
            LOG(fmt::format("Neumann Boundary {}, Function: {}",neumann_bc_tags[i],neumann_bc_s[i]));
        }
        LOG_VAR("Initial state function", init_state_s)
        LOG_VAR("Polynomial Degree", poly_deg)
        LOG_VAR("Finite Elements Type",fe_type)
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
