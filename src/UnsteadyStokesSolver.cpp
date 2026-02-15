#include <Generics.hpp>
#include <PDEParamsHandler.hpp>
#include <NavierStokesHelpers.hpp>
#include <UnsteadyStokesSolver.hpp>

namespace NavierStokesPDE{
    using namespace dealii;

    template<int dim>
    double UnsteadyStokesSolver<dim>::setup(String parameter_filename) {
        this->init(parameter_filename);

        LOG_TITLE("Reading Grid Serially")
        Triangulation<dim> non_parallel_mesh;
        GridIn<dim> grid_input;
        grid_input.attach_triangulation(non_parallel_mesh);

        std::ifstream mesh_file(this->mesh_filename);
        grid_input.read_msh(mesh_file);

        LOG_VAR("Cells",non_parallel_mesh.n_active_cells())

        LOG_TITLE("Distributing Grid to processors")
        GridTools::partition_triangulation(mpi_size, non_parallel_mesh);
        const auto construction_data =
            TriangulationDescription::Utilities::create_description_from_triangulation(
                non_parallel_mesh,
                MPI_COMM_WORLD
            );
        mesh.create_triangulation(construction_data);
        LOG_VAR("Cells distributed",mesh.n_global_active_cells())


        Timer timer;
        timer.start();


        LOG_TITLE("Init FE")
        LOG_VAR("Polynomial degree velocity",this->poly_deg_vel);
        LOG_VAR("Polynomial degree pressure",this->poly_deg_pressure);
        const FE_SimplexP<dim> fe_vel(this->poly_deg_vel);
        const FE_SimplexP<dim> fe_pressure(this->poly_deg_pressure);
        fe = std::make_unique<FESystem<dim>>(
            fe_vel, dim,
            fe_pressure, 1
        );
        LOG_VAR("DoF of FE System",fe->dofs_per_cell);


        LOG_TITLE("Init Quadrature")
        LOG_VAR("Quadrature points",fe->degree+1)
        quadrature = std::make_unique<QGaussSimplex<dim>>(fe->degree+1);
        boundary_quadrature = std::make_unique<QGaussSimplex<dim - 1>>(fe->degree+1);
        LOG_VAR("Quadrature points/cells",quadrature->size());
        LOG_VAR("Quadrature points/cells on boundaries",boundary_quadrature->size());


        LOG_TITLE("Init DoF")
        dof_handler.reinit(mesh);
        dof_handler.distribute_dofs(*fe);

        // have the components be [ux,uy,uz,p]
        std::vector<unsigned int> block_component(dim + 1, 0);
        block_component[dim] = 1;
        DoFRenumbering::component_wise(dof_handler, block_component);

        locally_owned_dofs = dof_handler.locally_owned_dofs();
        locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);

        std::vector<types::global_dof_index> dofs_per_block =
            DoFTools::count_dofs_per_fe_block(dof_handler, block_component);

        const unsigned int n_u = dofs_per_block[0];
        const unsigned int n_p = dofs_per_block[1];

        block_owned_dofs.resize(2);
        block_relevant_dofs.resize(2);

        block_owned_dofs[0]    = locally_owned_dofs.get_view(0, n_u);
        block_owned_dofs[1]    = locally_owned_dofs.get_view(n_u, n_u + n_p);

        block_relevant_dofs[0] = locally_relevant_dofs.get_view(0, n_u);
        block_relevant_dofs[1] = locally_relevant_dofs.get_view(n_u, n_u + n_p);

        LOG_VAR("Total DoFs",(n_u+n_p))
        LOG_VAR("Velocity DoFs",n_u)
        LOG_VAR("Pressure DoFs",n_p)


        LOG_TITLE("Init Algebraic structure")
        LOG_TITLE("Creating Sparsity Pattern")

        // u interacts with v (∇u*∇v)
        // p interacts with v (p*div(v))
        // u interacts with q (q*div(v))
        // not p and q
        // we use coupling to signal this
        Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1);
        for (unsigned int c = 0; c < dim + 1; ++c) {
            for (unsigned int d = 0; d < dim + 1; ++d) {
                if (c == dim && d == dim) coupling[c][d] = DoFTools::none;
                else coupling[c][d] = DoFTools::always;
            }
        }

        TrilinosWrappers::BlockSparsityPattern stiff_sp(
            block_owned_dofs,
            MPI_COMM_WORLD
        );
        DoFTools::make_sparsity_pattern(dof_handler, coupling ,stiff_sp);
        stiff_sp.compress();

        // reinitialize to build mass matrix
        // now only p and q
        for (unsigned int c = 0; c < dim + 1; ++c) {
            for (unsigned int d = 0; d < dim + 1; ++d) {
                if (c == dim && d == dim) coupling[c][d] = DoFTools::always;
                else coupling[c][d] = DoFTools::none;
            }
        }

        TrilinosWrappers::BlockSparsityPattern mass_sp(
            block_owned_dofs,
            MPI_COMM_WORLD
        );
        DoFTools::make_sparsity_pattern(dof_handler, coupling ,mass_sp);
        mass_sp.compress();


        LOG_TITLE("Creating Structure from Sparsity Pattern")
        stiff.reinit(stiff_sp);
        mass_pressure.reinit(mass_sp);

        rhs.reinit(block_owned_dofs, MPI_COMM_WORLD);
        owned_sol.reinit(block_owned_dofs, MPI_COMM_WORLD);
        sol.reinit(block_owned_dofs,block_relevant_dofs,MPI_COMM_WORLD);

        timer.stop();
        return timer.wall_time();
    }

    template<int dim>
    double UnsteadyStokesSolver<dim>::assemble() {
        Timer timer;
        timer.start();

        LOG_TITLE("Init Assembly")
        const auto & dofs_per_cell = fe->dofs_per_cell;
        const auto & n_q = quadrature->size();
        const unsigned int n_q_face      = boundary_quadrature->size();

        LOG_TITLE("Init Fe and Quadrature evaluator")

        FEValues<dim> fe_v(
            *fe,
            *quadrature,
            update_values |                 // \phi(x)
            update_gradients |              // \nabla\phi(x)
            update_quadrature_points |      // x_j \forall 1\leq j \leq n_q
            update_JxW_values               // J(x_j) and w_j
        );
        FEFaceValues<dim> fe_f_v(
            *fe,
            *boundary_quadrature,
            update_values |                 // \phi(x)
            update_normal_vectors |         // n at the border
            update_quadrature_points |      // x_j \forall 1\leq j \leq n_q
            update_JxW_values               // J(x_j) and w_j
        );

        LOG_TITLE("Init Local Algebraic System")
        FullMatrix<double> stiff_local(
            dofs_per_cell,
            dofs_per_cell
            );
        FullMatrix<double> mass_local(
            dofs_per_cell,
            dofs_per_cell
        );
        Vector<double> rhs_local(dofs_per_cell);

        // map to global indices
        // dof_indices(local) == global
        std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

        // clear global in case we called the function again;
        stiff = 0.0;
        mass_pressure = 0.0;
        rhs = 0.0;

        FEValuesExtractors::Vector vel(0); // 0 is the first component of the vector
        FEValuesExtractors::Scalar press(dim);

        std::vector<Tensor<1, dim>> current_velocity_values(n_q);
        std::vector<Tensor<2, dim>> current_velocity_gradients(n_q);
        std::vector<double> current_velocity_divergences(n_q);
        std::vector<double> current_pressure_values(n_q);

        LOG_TITLE("Filling Local Algebraic System")
        for (const auto & cell : dof_handler.active_cell_iterators()) {
            if (!cell->is_locally_owned()) continue;

            // computing requested values (flags) on cell
            fe_v.reinit(cell);

            // clear old matrix and rhs to compute on new cell
            stiff_local = 0.0;
            rhs_local   = 0.0;
            mass_local  = 0.0;

            fe_v[vel].get_function_values(sol,current_velocity_values);
            fe_v[vel].get_function_gradients(sol, current_velocity_gradients);
            fe_v[vel].get_function_divergences(sol, current_velocity_divergences);
            fe_v[press].get_function_values(sol,current_pressure_values);

            for (unsigned int q=0; q<n_q; ++q) {
                for (unsigned int i=0; i<dofs_per_cell; ++i) {
                    for (unsigned int j=0; j<dofs_per_cell; ++j) {
                        // M/Δt
                        stiff_local(i,j) += (1/this->dt) * (
                                scalar_product(
                                fe_v[vel].value(i,q),
                                fe_v[vel].value(j,q)
                                )
                            ) * fe_v.JxW(q);
                        // + θ*a(u_h,v_h)
                        stiff_local(i,j) += this->theta * (
                                (this->nu * scalar_product(fe_v[vel].gradient(i,q),fe_v[vel].gradient(j,q))) +
                                (this->sigma * scalar_product(fe_v[vel].value(i,q),fe_v[vel].value(j,q)))
                            ) * fe_v.JxW(q);
                        // + θ*b(v,p) = (-∫ p⋅div(v) dx)
                        stiff_local(i,j) -= this->theta *
                                            fe_v[vel].divergence(i,q) *
                                            fe_v[press].value(j,q) *
                                            fe_v.JxW(q);
                        // + θ*b(u,q)
                        stiff_local(i,j) -= this->theta *
                                            fe_v[vel].divergence(j,q) *
                                            fe_v[press].value(i,q) *
                                            fe_v.JxW(q);

                        // M_ij = psy_i*psy_j
                        mass_local(i,j) +=  fe_v[press].value(i,q) *
                                            fe_v[press].value(j,q) /
                                            this->nu * fe_v.JxW(q) ;
                    }

                    // M/Δt
                    rhs_local(i) += (1/this->dt) * (
                            scalar_product(
                            fe_v[vel].value(i,q),
                            current_velocity_values[q]
                            )
                        ) * fe_v.JxW(q);

                    // + (1-θ)*a(u_h,v_h)
                    rhs_local(i) -= (1-this->theta) * (
                            (this->nu * scalar_product(fe_v[vel].gradient(i,q),current_velocity_gradients[q])) +
                            (this->sigma * scalar_product(fe_v[vel].value(i,q),current_velocity_values[q]))
                        ) * fe_v.JxW(q);
                    // + θ*b(v,p) = (-∫ p⋅div(v) dx)
                    rhs_local(i) += (1-this->theta) *
                                        fe_v[vel].divergence(i,q) *
                                        current_pressure_values[q] *
                                        fe_v.JxW(q);
                    // + θ*b(u,q)
                    rhs_local(i) += (1-this->theta) *
                                        current_velocity_divergences[q] *
                                        fe_v[press].value(i,q) *
                                        fe_v.JxW(q);

                    // this remains as is, if f(x,t) = f(x) then: θF + (1-θ)F = F
                    Tensor<1,dim> force_tensor;
                    for (auto d=0; d<dim; ++d) {
                        force_tensor[d] = force.value(fe_v.quadrature_point(q),d);
                    }
                    rhs_local(i) += force_tensor *
                                    fe_v[vel].value(i,q) *
                                    fe_v.JxW(q);
                }
            }

            if (cell->at_boundary()) {
                for (unsigned int face_id = 0 ; face_id < cell->n_faces(); ++face_id) {
                    auto face = cell->face(face_id);
                    auto tag = face->boundary_id();
                    if (face->at_boundary() && (tag == 0 || tag == 1 /* && other tags*/)) {
                        fe_f_v.reinit(cell,face_id);
                        for (unsigned int q=0; q<n_q_face; ++q) {
                            for (unsigned int i = 0; i<dofs_per_cell; ++i) {
                                // this remains as is
                                // if ψ(x,t) = ψ(x) then: θψ + (1-θ)ψ = ψ
                                rhs_local(i) +=
                                    (tag == 0 ? n_bc_0.scalar_value() : n_bc_1.scalar_value()) *
                                    scalar_product(
                                        fe_f_v.normal_vector(q),
                                        fe_f_v[vel].value(i,q)
                                        ) *
                                    fe_f_v.JxW(q);
                            }
                        }
                    }
                }
            }

            // retrieve local to global indices conversion for this cell
            cell->get_dof_indices(dof_indices);

            // reduce
            stiff.add(dof_indices, stiff_local);
            rhs.add(dof_indices, rhs_local);
            mass_pressure.add(dof_indices, mass_local);
        }

        LOG_TITLE("Synchronizing computation of Algebraic systems")
        stiff.compress(VectorOperation::add);
        rhs.compress(VectorOperation::add);
        mass_pressure.compress(VectorOperation::add);

        LOG_TITLE("Applying Boundary Conditions")
        // boundary_vals(loc_index) == phi(loc_index)
        BoundaryMap boundary_vals;
        // associate boundaryID to function
        BoundaryFunctionMap boundary_functions;
        // Boundary function will be defined for 4 dimensions, but we mask the last one
        ComponentMask mask_velocity(dim + 1, true);
        // turn of dimension 4 (3 starting from 0)
        mask_velocity.set(dim, false);

        // boundary_functions[0] = &d_bc_0;
        // // boundary_functions[1] = &d_bc_1;
        // // boundary_functions[2] = &d_bc_2;
        // // ...
        //
        // // populate map boundary_vals
        // VectorTools::interpolate_boundary_values(
        //     dof_handler,
        //     boundary_functions,
        //     boundary_vals,
        //     mask_velocity
        // );

        LOG_TITLE("Applying Wall Conditions")
        // wall conditions (aka homogeneous dirichlet) should go last
        boundary_functions.clear();
        Functions::ZeroFunction<dim> zero(dim+1);
        boundary_functions[2] = &zero;
        boundary_functions[3] = &zero;

        // populate map boundary_vals
        VectorTools::interpolate_boundary_values(
            dof_handler,
            boundary_functions,
            boundary_vals,
            mask_velocity
        );

        // modify matrix to respect boundary conditions
        MatrixTools::apply_boundary_values(
            boundary_vals,
            stiff,
            owned_sol,
            rhs,
            false
        );

        timer.stop();
        return timer.wall_time();
    }

    template<int dim>
    double UnsteadyStokesSolver<dim>::solve() {
        Timer timer;
        timer.start();

        LOG_TITLE("Initializing Solver")
        LOG_VAR("Max number of Iterations", this->max_iters)
        LOG_VAR("Tolerance", this->epsilon)
        // set all three values for the solver to use
        SolverControl solve_controller(
            this->max_iters,
            this->epsilon * rhs.l2_norm()
        );

        LOG_TITLE("Using BlockTriangular preconditioner")
        NavierStokesPDE::Preconditioners::BlockTriangular preconditioner;
        preconditioner.initialize(
            stiff.block(0,0),
            mass_pressure.block(1,1),
            stiff.block(1,0)
        );

        LOG_TITLE("Using GMRES solver")
        SolverGMRES<BlkVector> gmres(solve_controller);
        LOG_TITLE("Solving...")
        gmres.solve(stiff, owned_sol, rhs, preconditioner);

        LOG_VAR("Iterations performed", solve_controller.last_step())

        sol = owned_sol;

        timer.stop();
        return timer.wall_time();
    }

    template<int dim>
    void UnsteadyStokesSolver<dim>::output() const {

        LOG_TITLE("Preparing output structure")
        DataOut<dim> d_o;
        // a std::vec of 4 components first we initialize 3 for velocity
        InterpretationVector interpretation(
            dim,
            DataComponentInterpretation::component_is_part_of_vector
        );
        // initialize the last for pressure
        interpretation.push_back(DataComponentInterpretation::component_is_scalar);

        // same for the names of the components
        std::vector<String> names(dim, "velocity");
        names.push_back("pressure");

        d_o.add_data_vector(dof_handler, sol, names, interpretation);

        LOG_TITLE("Preparing mesh subdivision vector")
        std::vector<unsigned int> partition_tags(mesh.n_active_cells());
        GridTools::get_subdomain_association(mesh, partition_tags);

        const Vector<double> partitioning(partition_tags.begin(), partition_tags.end());
        d_o.add_data_vector(partitioning, "partitioning");

        LOG_TITLE("Finalizing output structure")
        d_o.build_patches();


        const std::filesystem::path mesh_path(this->output_filename);
        String output_dir = mesh_path.parent_path().string();
        String output_name = mesh_path.stem().string();
        LOG_VAR("Opening directory", output_dir)
        LOG_VAR("Opening Files Named", output_name)
        String out_file = d_o.write_vtu_with_pvtu_record(
            output_dir + "/",
            output_name,
            time_step,
            MPI_COMM_WORLD
        );
        LOG_TITLE("Successful Write")
        LOG_VAR("Logs written to",out_file)
    }

    template<int dim>
   double UnsteadyStokesSolver<dim>::run(const String &parameter_filename) {
        Timer timer;
        timer.start();

        const bool print_iter = true;

        auto set_up_time = setup(parameter_filename);
        LOG_ANY("Setup time: {:0.4f}s", set_up_time);

        const double T = this->T;
        const double dt = this->dt;
        const int expected_steps = static_cast<int>(T/dt);

        LOG_TITLE("Forcing Initial state")
        VectorTools::interpolate(
            dof_handler,
            initial_sol,
            owned_sol
        );
        sol = owned_sol;    // share with processors

        time = 0.0;
        time_step = 0;

        output();

        LOG_TITLE("Starting Time Loop")
        LOG_VAR("Expected number of cycles", expected_steps)

        while (time < (T - 0.5*dt)) {
            time += dt;
            ++time_step;
            if (print_iter)
                pcout << fmt::format("{:<90}",fmt::format("{} == {}/{}","Starting Time step", time_step,expected_steps)) << std::endl;
            double assemble_time = assemble();
            LOG_ANY("RHS assemble time: {:0.4f}s", assemble_time);
            double solve_time = solve();
            LOG_ANY("Solve time: {:0.4f}s", solve_time);

            output();
        }

        timer.stop();
        return timer.wall_time();
    }

    /*
    template<int dim>
    double UnsteadyStokesSolver<dim>::compare_solution(
        const VectorTools::NormType &norm_type,
        const Function<dim> &exact
    ) const {
        // init quadrature to intergrate error
        const QGaussSimplex<dim> quadrature_error(this->poly_deg + 2);
        // triangular mesh needs mapping
        FE_SimplexP<dim> linear_fe(1);
        MappingFE        mapping(linear_fe);
        // init error vector on the cells
        Vector<double> error_on_cell(mesh.n_active_cells());
        // compute error
        VectorTools::integrate_difference(
            mapping,
            dof_handler,
            sol,
            exact,
            error_on_cell,
            quadrature_error,
            norm_type
        );
        return VectorTools::compute_global_error(mesh,error_on_cell,norm_type);
    }
    */

}
