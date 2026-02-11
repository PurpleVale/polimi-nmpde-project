#include <Generics.hpp>
#include <PDEParamsHandler.hpp>

namespace NavierStokesPDE {
    using namespace dealii;

    template<int dim>
    void StokesParamHandler<dim>::declare_parameters() {
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
        prm.enter_subsection("Finite Elements");
        {
            prm.declare_entry(
                "Polynomial degree velocity",
                "2",
                Patterns::Integer(1,7),
                "The polynomial degree of the Elements of the velocity"
            );
            prm.declare_entry(
                "Polynomial degree pressure",
                "1",
                Patterns::Integer(1,7),
                "The polynomial degree of the Elements of the pressure"
            );
        }
        prm.leave_subsection();
        prm.enter_subsection("Coefficients");
        {
            prm.declare_entry(
                "Kinematic Viscosity",
                "1.0",
                Patterns::Double(1e-16),
                "The Kinematic Viscosity Coefficient (nu)"
            );
            prm.declare_entry(
                "Sigma",
                "0.0",
                Patterns::Double(0.0),
                "The coefficient sigma"
            );
        }
        prm.leave_subsection();
        param_initialized = true;
    }

    template<int dim>
    void StokesParamHandler<dim>::init(const String &filename) {
        if (initialized) return;

        if (!param_initialized) {
            declare_parameters();
        }

        LOG_TITLE("Initializing Elliptic parameters");
        constants["pi"] = numbers::PI;

        LOG_TITLE("Reading Parameter File")
        LOG_VAR("File Name", filename)
        prm.parse_input(filename);

        prm.enter_subsection("Solver");
        {
            LOG_TITLE("Getting Solver Parameters")
            max_iters = prm.get_integer("Max iters");
            epsilon = prm.get_double("Tolerance");
        }
        prm.leave_subsection();
        prm.enter_subsection("Files");
        {
            LOG_TITLE("Getting File Names")
            mesh_filename = prm.get("Mesh file");
            output_filename = prm.get("Output file");
        }
        prm.leave_subsection();
        prm.enter_subsection("Finite Elements");
        {
            LOG_TITLE("Getting Finite Elements Parameters")
            poly_deg_vel = prm.get_integer("Polynomial degree velocity");
            poly_deg_pressure = prm.get_integer("Polynomial degree pressure");
        }
        prm.leave_subsection();
        prm.enter_subsection("Coefficients");
        {
            LOG_TITLE("Getting Coefficients")
            nu = prm.get_double("Kinematic Viscosity");
            sigma = prm.get_double("Sigma");
        }
        prm.leave_subsection();

        LOG_TITLE("Parameters Read successful")
        LOG_VAR("Max iters",max_iters)
        LOG_VAR("Tolerance",epsilon)
        LOG_VAR("Mesh file",mesh_filename)
        LOG_VAR("Output file",output_filename)
        LOG_VAR("Polynomial Degree for velocity", poly_deg_vel)
        LOG_VAR("Polynomial Degree for pressure", poly_deg_pressure)
        LOG_VAR("Kinematic Viscosity",nu);
        LOG_VAR("Sigma",sigma);
        initialized = true;
    }

    template<int dim>
    void StokesParamHandler<dim>::print_parameters(const String &filename) {
        print_parameters_as(filename,ParameterHandler::DefaultStyle);
    }

    template<int dim>
    void StokesParamHandler<dim>::print_editable_parameters(const String &filename) {
        print_parameters_as(filename,ParameterHandler::XML);
    }

    template<int dim>
    void StokesParamHandler<dim>::print_parameters_as(
        const String &filename,
        ParameterHandler::OutputStyle style
    ) {
        if (! param_initialized) {
            declare_parameters();
        }
        prm.print_parameters(filename,style);
    }


    template<int dim>
    void UnsteadyStokesParamHandler<dim>::declare_parameters() {
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
        prm.enter_subsection("Finite Elements");
        {
            prm.declare_entry(
                "Polynomial degree velocity",
                "2",
                Patterns::Integer(1,7),
                "The polynomial degree of the Elements of the velocity"
            );
            prm.declare_entry(
                "Polynomial degree pressure",
                "1",
                Patterns::Integer(1,7),
                "The polynomial degree of the Elements of the pressure"
            );
        }
        prm.leave_subsection();
        prm.enter_subsection("Coefficients");
        {
            prm.declare_entry(
                "Kinematic Viscosity",
                "1.0",
                Patterns::Double(1e-16),
                "The Kinematic Viscosity Coefficient (nu)"
            );
            prm.declare_entry(
                "Sigma",
                "0.0",
                Patterns::Double(0.0),
                "The coefficient sigma"
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
    void UnsteadyStokesParamHandler<dim>::init(const String &filename) {
        if (initialized) return;

        if (!param_initialized) {
            declare_parameters();
        }

        LOG_TITLE("Initializing Elliptic parameters");
        constants["pi"] = numbers::PI;

        LOG_TITLE("Reading Parameter File")
        LOG_VAR("File Name", filename)
        prm.parse_input(filename);

        prm.enter_subsection("Solver");
        {
            LOG_TITLE("Getting Solver Parameters")
            max_iters = prm.get_integer("Max iters");
            epsilon = prm.get_double("Tolerance");
        }
        prm.leave_subsection();
        prm.enter_subsection("Files");
        {
            LOG_TITLE("Getting File Names")
            mesh_filename = prm.get("Mesh file");
            output_filename = prm.get("Output file");
        }
        prm.leave_subsection();
        prm.enter_subsection("Finite Elements");
        {
            LOG_TITLE("Getting Finite Elements Parameters")
            poly_deg_vel = prm.get_integer("Polynomial degree velocity");
            poly_deg_pressure = prm.get_integer("Polynomial degree pressure");
        }
        prm.leave_subsection();
        prm.enter_subsection("Coefficients");
        {
            LOG_TITLE("Getting Coefficients")
            nu = prm.get_double("Kinematic Viscosity");
            sigma = prm.get_double("Sigma");
        }
        prm.leave_subsection();
        prm.enter_subsection("Time Parameters");
        {
            LOG_TITLE("Getting Time Parameters")
            T = prm.get_double("Final time");
            dt = prm.get_double("Time step");
            theta = prm.get_double("Theta");
        }
        prm.leave_subsection();

        LOG_TITLE("Parameters Read successful")
        LOG_VAR("Max iters",max_iters)
        LOG_VAR("Tolerance",epsilon)
        LOG_VAR("Mesh file",mesh_filename)
        LOG_VAR("Output file",output_filename)
        LOG_VAR("Polynomial Degree for velocity", poly_deg_vel)
        LOG_VAR("Polynomial Degree for pressure", poly_deg_pressure)
        LOG_VAR("Kinematic Viscosity",nu);
        LOG_VAR("Sigma",sigma);
        LOG_VAR("Time limit of the PDE",T)
        LOG_VAR("Time step",dt)
        LOG_VAR("Theta Method Parameter",theta)
        initialized = true;
    }

    template<int dim>
    void UnsteadyStokesParamHandler<dim>::print_parameters(const String &filename) {
        print_parameters_as(filename,ParameterHandler::DefaultStyle);
    }

    template<int dim>
    void UnsteadyStokesParamHandler<dim>::print_editable_parameters(const String &filename) {
        print_parameters_as(filename,ParameterHandler::XML);
    }

    template<int dim>
    void UnsteadyStokesParamHandler<dim>::print_parameters_as(
        const String &filename,
        ParameterHandler::OutputStyle style
    ) {
        if (! param_initialized) {
            declare_parameters();
        }
        prm.print_parameters(filename,style);
    }

}
