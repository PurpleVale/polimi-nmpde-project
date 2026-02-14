#include <Generics.hpp>
#include <PDEParamsHandler.hpp>

namespace EllipticPDE {

    template<int dim>
    void DNParamHandler<dim>::declare_parameters() {
        LOG_TITLE("Constructing Parameters Structure")
        prm.enter_subsection("Iteration controls");
        {
            prm.declare_entry(
                "Max iters",
                "100",
                Patterns::Integer(1,100000),
                "Max number of iters"
            );
            prm.declare_entry(
                "Tolerance",
                "1e-4",
                Patterns::Double(),
                "Tolerance for the increment of the Neumann subdomain"
            );
        }
        prm.leave_subsection();
        prm.enter_subsection("Files");
        {
            prm.declare_entry(
                "Dirichlet Parameter file",
                "./dirichlet.xml",
                Patterns::FileName(Patterns::FileName::input),
                "A .xml file containing the parameters for the first subdomain"
            );
            prm.declare_entry(
                "Neumann Parameter file",
                "./neumann.xml",
                Patterns::FileName(Patterns::FileName::input),
                "A .xml file containing the parameters for the second subdomain"
            );
        }
        prm.leave_subsection();
        prm.declare_entry(
            "Relaxation",
            "1.0",
            Patterns::Double(),
            "The relaxation parameter, 1 is no relaxation"
        );
        param_initialized = true;
    }

    template<int dim>
    void DNParamHandler<dim>::init(const String &filename) {
        if (initialized) return;

        if (!param_initialized) {
            declare_parameters();
        }

        LOG_TITLE("Initializing DN parameters");

        LOG_TITLE("Reading Parameter File")
        LOG_VAR("File Name", filename)
        prm.parse_input(filename);

        prm.enter_subsection("Iteration controls");
        {
            max_iters = prm.get_integer("Max iters");
            tolerance = prm.get_double("Tolerance");
        }
        prm.leave_subsection();
        prm.enter_subsection("Files");
        {
            dirichlet_prm_file = prm.get("Dirichlet Parameter file");
            neumann_prm_file = prm.get("Neumann Parameter file");
        }
        prm.leave_subsection();

        relaxation = prm.get_double("Relaxation");

        LOG_TITLE("DN Parameters Read successful")
        LOG_VAR("Max number of iterations",max_iters)
        LOG_VAR("Tolerance",tolerance)
        LOG_VAR("Dirichlet Parameter file",dirichlet_prm_file)
        LOG_VAR("Neumann Parameter file",neumann_prm_file)
        LOG_VAR("Relaxation",relaxation)
        initialized = true;
    }

    template<int dim>
    void DNParamHandler<dim>::print_parameters(const String &filename) {
        print_parameters_as(filename,ParameterHandler::DefaultStyle);
    }

    template<int dim>
    void DNParamHandler<dim>::print_editable_parameters(const String &filename) {
        print_parameters_as(filename,ParameterHandler::XML);
    }

    template<int dim>
    void DNParamHandler<dim>::print_parameters_as(
        const String &filename,
        ParameterHandler::OutputStyle style
    ) {
        if (! param_initialized) {
            declare_parameters();
        }
        prm.print_parameters(filename,style);
    }
}
