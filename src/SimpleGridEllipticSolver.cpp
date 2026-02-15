#include <Generics.hpp>
#include <PDEParamsHandler.hpp>
#include <SimpleGridEllipticSolver.hpp>

namespace EllipticPDE{
    using namespace dealii;

    double SimpleGridEllipticSolver::setup(String parameter_filename) {
        this->init(parameter_filename);

        LOG_TITLE("Creating Grid")
        Triangulation<dim> non_parallel_mesh;
        GridGenerator::subdivided_hyper_cube(non_parallel_mesh, N_elm, 0.0, 1.0, true);

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
        LOG_VAR("FE type",this->fe_type);
        LOG_VAR("Polynomial degree",this->poly_deg);
        if (this->fe_type == "Q") {
            fe = std::make_unique<FE_Q<dim>>(this->poly_deg);
        } else {
            fe = std::make_unique<FE_SimplexP<dim>>(this->poly_deg);
        }
        LOG_VAR("DoF of FE",fe->dofs_per_cell);

        LOG_TITLE("Init Quadrature")
        LOG_VAR("Quadrature points",this->poly_deg + 1)
        if (this->fe_type == "Q") {
            quadrature = std::make_unique<QGauss<dim>>(this->poly_deg+1);
            boundary_quadrature = std::make_unique<QGauss<dim - 1>>(this->poly_deg + 1);
        } else {
            quadrature = std::make_unique<QGaussSimplex<dim>>(this->poly_deg+1);
            boundary_quadrature = std::make_unique<QGaussSimplex<dim - 1>>(this->poly_deg + 1);
        }
        LOG_VAR("Quadrature points/cells",quadrature->size());
        LOG_VAR("Quadrature points/cells on boundaries",boundary_quadrature->size());

        LOG_TITLE("Init DoF")
        dof_handler.reinit(mesh);
        dof_handler.distribute_dofs(*fe);
        locally_owned_dofs = dof_handler.locally_owned_dofs();
        LOG_VAR("number of DoFs",dof_handler.n_dofs())

        LOG_TITLE("Init Algebraic structure")
        LOG_TITLE("Creating Sparsity Pattern")
        TrilinosWrappers::SparsityPattern stiff_sp(
            locally_owned_dofs,
            MPI_COMM_WORLD
        );
        DoFTools::make_sparsity_pattern(dof_handler, stiff_sp);
        stiff_sp.compress();

        LOG_TITLE("Creating Structure from Sparsity Pattern")
        stiff.reinit(stiff_sp);

        rhs.reinit(locally_owned_dofs, MPI_COMM_WORLD);
        sol.reinit(locally_owned_dofs, MPI_COMM_WORLD);

        timer.stop();
        return timer.wall_time();
    }

    double SimpleGridEllipticSolver::assemble() {
        Timer timer;
        timer.start();

        LOG_TITLE("Init Assembly")
        const auto & dofs_per_cell = fe->dofs_per_cell;
        const auto & n_q = quadrature->size();
        LOG_VAR("Dof/cell",dofs_per_cell)
        LOG_VAR("Quadrature points", n_q)

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
            update_quadrature_points |      // x_j \forall 1\leq j \leq n_q
            update_JxW_values               // J(x_j) and w_j
        );

        LOG_TITLE("Init Local Algebraic System")
        FullMatrix<double> stiff_local(
            dofs_per_cell,
            dofs_per_cell
            );
        stiff_local = 0.0; // make sure it's empty
        Vector<double> rhs_local(dofs_per_cell);
        rhs_local = 0.0; // make sure it's empty

        // map to global indices
        // dof_indices(local) == global
        std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

        // clear global in case we called the function again;
        stiff = 0.0;
        rhs = 0.0;

        LOG_TITLE("Filling Local Algebraic System")
        for (const auto & cell : dof_handler.active_cell_iterators()) {
            if (!cell->is_locally_owned()) continue;

            // clear old matrix and rhs to compute on new cell
            stiff_local = 0.0;
            rhs_local = 0.0;

            // computing requested values (flags) on cell
            fe_v.reinit(cell);

            for (unsigned int q=0; q<n_q; ++q) {
                auto x_q = fe_v.quadrature_point(q);
                // iterate over \phi_i/\phi_j on the local cell
                for (unsigned int i=0; i<dofs_per_cell; ++i) {
                    for (unsigned int j=0; j<dofs_per_cell; ++j) {
                        Tensor<1,dim> advection_tensor;
                        for (auto d = 0; d<dim ; ++d)
                            advection_tensor[d] = this->advection_c.value(x_q,d);
                        stiff_local(i,j) += (
                                (
                                    this->diffusion_c.value(x_q) *
                                    fe_v.shape_grad(i,q) *
                                    fe_v.shape_grad(j,q)
                                ) + (
                                    advection_tensor *
                                    fe_v.shape_grad(j,q) *
                                    fe_v.shape_value(i,q)
                                ) + (
                                    this->reaction_c.value(x_q) *
                                    fe_v.shape_value(i,q) *
                                    fe_v.shape_value(j,q)
                                )
                            ) *  fe_v.JxW(q);
                    }
                rhs_local(i) += this->force_term.value(x_q) * fe_v.shape_value(i,q) * fe_v.JxW(q);
                }
            }

            if (cell->at_boundary()) {
                for (unsigned int face_id = 0 ; face_id < cell->n_faces(); ++face_id) {
                    auto face = cell->face(face_id);
                    auto pos = std::find(this->neumann_bc_tags.begin(),this->neumann_bc_tags.end(),face->boundary_id());
                    if (face->at_boundary() && pos != this->neumann_bc_tags.end()) {
                        fe_f_v.reinit(cell,face_id);
                        auto neumann_index = pos - this->neumann_bc_tags.begin();
                        for (unsigned int q=0; q<boundary_quadrature->size(); ++q) {
                            auto x_q = fe_f_v.quadrature_point(q);
                            for (unsigned int i = 0; i<dofs_per_cell; ++i) {
                                rhs_local(i) += this->neumann_bc[neumann_index]->value(x_q) * fe_f_v.shape_value(i,q) * fe_f_v.JxW(q);
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
        }

        LOG_TITLE("Synchronizing computation of Algebraic systems")
        stiff.compress(VectorOperation::add);
        rhs.compress(VectorOperation::add);

        LOG_TITLE("Applying Boundary Conditions")
        // boundary_vals(loc_index) == phi(loc_index)
        BoundaryMap boundary_vals;
        // associate boundaryID to function
        BoundaryFunctionMap boundary_functions;
        for (long unsigned int i = 0; i < this->dirichlet_bc_tags.size(); i++)
            boundary_functions[this->dirichlet_bc_tags[i]] = this->dirichlet_bc[i].get();

        // populate map boundary_vals
        VectorTools::interpolate_boundary_values(
            dof_handler,
            boundary_functions,
            boundary_vals
        );

        // modify matrix to respect boundary conditions
        MatrixTools::apply_boundary_values(
            boundary_vals,
            stiff,
            sol,
            rhs,
            false
        );

        timer.stop();
        return timer.wall_time();
    }

    double SimpleGridEllipticSolver::solve() {
        Timer timer;
        timer.start();

        LOG_TITLE("Initializing Solver")
        LOG_VAR("Max number of Iterations", this->max_iters)
        LOG_VAR("Absolute precision", this->epsilon)
        LOG_VAR("Relative precision (residual)", this->residual_tolerance)
        // set all three values for the solver to use
        ReductionControl solve_controller(
            this->max_iters,
            this->epsilon,
            this->residual_tolerance
        );

        if (this->preconditioner == "Jacobi") {
            LOG_TITLE("Using Jacobi preconditioner")
            TrilinosWrappers::PreconditionJacobi prec;
            prec.initialize(stiff,TrilinosWrappers::PreconditionJacobi::AdditionalData(0.6));
            if (this->symmetric_solver) {
                LOG_TITLE("Using CG solver")
                SolverCG<TrilinosWrappers::MPI::Vector> cg(solve_controller);
                LOG_TITLE("Solving...")
                cg.solve(stiff, sol, rhs, prec);
            } else {
                LOG_TITLE("Using GMRES solver")
                SolverGMRES<TrilinosWrappers::MPI::Vector> gmres(solve_controller);
                LOG_TITLE("Solving...")
                gmres.solve(stiff, sol, rhs, prec);
            }
        } else if (this->preconditioner == "SSOR") {
            LOG_TITLE("Using SSOR preconditioner")
            TrilinosWrappers::PreconditionSSOR prec;
            prec.initialize(stiff,TrilinosWrappers::PreconditionSSOR::AdditionalData(1.0));
            if (this->symmetric_solver) {
                LOG_TITLE("Using CG solver")
                SolverCG<TrilinosWrappers::MPI::Vector> cg(solve_controller);
                LOG_TITLE("Solving...")
                cg.solve(stiff, sol, rhs, prec);
            } else {
                LOG_TITLE("Using GMRES solver")
                SolverGMRES<TrilinosWrappers::MPI::Vector> gmres(solve_controller);
                LOG_TITLE("Solving...")
                gmres.solve(stiff, sol, rhs, prec);
            }
        } else if (this->preconditioner == "ILU") {
            LOG_TITLE("Using ILU preconditioner")
            TrilinosWrappers::PreconditionILU prec;
            prec.initialize(stiff,TrilinosWrappers::PreconditionILU::AdditionalData(2));
            if (this->symmetric_solver) {
                LOG_TITLE("Using CG solver")
                SolverCG<TrilinosWrappers::MPI::Vector> cg(solve_controller);
                LOG_TITLE("Solving...")
                cg.solve(stiff, sol, rhs, prec);
            } else {
                LOG_TITLE("Using GMRES solver")
                SolverGMRES<TrilinosWrappers::MPI::Vector> gmres(solve_controller);
                LOG_TITLE("Solving...")
                gmres.solve(stiff, sol, rhs, prec);
            }
        } else if (this->preconditioner == "ILUT") {
            LOG_TITLE("Using ILUT preconditioner")
            TrilinosWrappers::PreconditionILUT prec;
            prec.initialize(stiff,TrilinosWrappers::PreconditionILUT::AdditionalData(0,2));
            if (this->symmetric_solver) {
                LOG_TITLE("Using CG solver")
                SolverCG<TrilinosWrappers::MPI::Vector> cg(solve_controller);
                LOG_TITLE("Solving...")
                cg.solve(stiff, sol, rhs, prec);
            } else {
                LOG_TITLE("Using GMRES solver")
                SolverGMRES<TrilinosWrappers::MPI::Vector> gmres(solve_controller);
                LOG_TITLE("Solving...")
                gmres.solve(stiff, sol, rhs, prec);
            }
        } else {
            LOG_TITLE("Not using a preconditioner")
            if (this->symmetric_solver) {
                LOG_TITLE("Using CG solver")
                SolverCG<TrilinosWrappers::MPI::Vector> cg(solve_controller);
                LOG_TITLE("Solving...")
                cg.solve(stiff, sol, rhs, TrilinosWrappers::PreconditionIdentity());
            } else {
                LOG_TITLE("Using GMRES solver")
                SolverGMRES<TrilinosWrappers::MPI::Vector> gmres(solve_controller);
                LOG_TITLE("Solving...")
                gmres.solve(stiff, sol, rhs, TrilinosWrappers::PreconditionIdentity());
            }
        }
        LOG_VAR("Iterations performed", solve_controller.last_step())

        timer.stop();
        return timer.wall_time();
    }

    void SimpleGridEllipticSolver::output() const {

        LOG_TITLE("Synchronizing with other processors")
        const IndexSet locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);
        // ghost are part of the solution computed that do not belong to me (borders)
        TrilinosWrappers::MPI::Vector sol_ghost(
            locally_owned_dofs,
            locally_relevant_dofs,
            MPI_COMM_WORLD
        );
        sol_ghost = sol; // perform comunication implicitly

        LOG_TITLE("Preparing output structure")
        DataOut<dim> d_o;
        // add ghost solution
        d_o.add_data_vector(dof_handler, sol_ghost, "solution");

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
            0,
            MPI_COMM_WORLD
        );
        LOG_TITLE("Successful Write")
        LOG_VAR("Logs written to",out_file)
    }

    double SimpleGridEllipticSolver::compare_solution(
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

}
