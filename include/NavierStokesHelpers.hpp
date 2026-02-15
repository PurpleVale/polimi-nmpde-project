#ifndef NNPDE_PROJECT_NAVIERSTOKESHELPERS_HPP
#define NNPDE_PROJECT_NAVIERSTOKESHELPERS_HPP

#include <Generics.hpp>
#include <deal.II/lac/trilinos_parallel_block_vector.h>

namespace NavierStokesPDE {
    using namespace dealii;

    using String = std::string;
    using BoundaryMap =  std::map<types::global_dof_index, double>;
    using BoundaryIds  = types::boundary_id;
    using BlkVector = TrilinosWrappers::MPI::BlockVector;
    using SpMtx = TrilinosWrappers::SparseMatrix;
    using ILU = TrilinosWrappers::PreconditionILU;
    using DVector = TrilinosWrappers::MPI::Vector;

    namespace StokesFunctions {

        template<int dim> const Tensor<1,dim,double> constant_force; // initialize else where

        template<int dim>
        class Force : public Function<dim> {
        public:
            Force(const unsigned int components) : Function<dim>(components) {}

            virtual
            void vector_value(const Point<dim> &p, Vector<double> &values) const override {
                values[0] = x_expr(p);
                if constexpr (dim >= 2) values[1] = y_expr(p);
                if constexpr (dim >= 3) values[2] = z_expr(p);
            }

            virtual
            double value(const Point<dim> &p, const unsigned int component = 0) const override {
                if (component == 0) return x_expr(p);
                if (component == 1) return y_expr(p);
                return z_expr(p);
            }

        protected:
            double x_expr(const Point<dim> &p) const {
                return 1.0 + std::sin(p[0])*std::sin(p[0]);
            }
            double y_expr(const Point<dim> &p) const {
                return 0.0;
            }
            double z_expr(const Point<dim> &p) const {
                return 0.0;
            }
        };

        template<int dim>
        class InitialSolution : public Function<dim> {
        public:
            InitialSolution(const unsigned int components) : Function<dim>(components) {}

            virtual
            void vector_value(const Point<dim> &p, Vector<double> &values) const override {
                values[0] = x_expr(p);
                if constexpr (dim >= 2) values[1] = y_expr(p);
                if constexpr (dim >= 3) values[2] = z_expr(p);
            }

            virtual
            double value(const Point<dim> &p, const unsigned int component = 0) const override {
                if (component == 0) return x_expr(p);
                if (component == 1) return y_expr(p);
                return z_expr(p);
            }

        protected:
            double x_expr(const Point<dim> &p) const {
                return 0.0;
            }
            double y_expr(const Point<dim> &p) const {
                return 0.0;
            }
            double z_expr(const Point<dim> &p) const {
                return 0.0;
            }
        };

        template<int dim>
        class DirichletCondition0 : public Function<dim> {
            public:
            DirichletCondition0() : Function<dim>(dim+1) {}

            virtual
            void vector_value(const Point<dim> &p, Vector<double> &values) const override {
                values[0] = x_expr(p);
                if constexpr (dim >= 2) values[1] = y_expr(p);
                if constexpr (dim >= 3) values[2] = z_expr(p);
            }

            virtual
            double value(const Point<dim> &p, const unsigned int component = 0) const override {
                if (component == 0) return x_expr(p);
                if (component == 1) return y_expr(p);
                return z_expr(p);
            }

            protected:
            double x_expr(const Point<dim> &p) const {
                return 0.0;
            }
            double y_expr(const Point<dim> &p) const {
                return 0.0;
            }
            double z_expr(const Point<dim> &p) const {
                return 0.0;
            }
        };

        template<int dim>
        class DirichletCondition1 : public Function<dim> {
            public:
                DirichletCondition1() : Function<dim>(dim + 1) {}

                virtual
                void vector_value(const Point<dim> &p, Vector<double> &values) const override {
                    values[0] = x_expr(p);
                    if constexpr (dim >= 2) values[1] = y_expr(p);
                    if constexpr (dim >= 3) values[2] = z_expr(p);
                }

                virtual
                double value(const Point<dim> &p, const unsigned int component = 0) const override {
                    if (component == 0) return x_expr(p);
                    if (component == 1) return y_expr(p);
                    return z_expr(p);
                }

                protected:
                double x_expr(const Point<dim> &p) const {
                    return 0.0;
                }
                double y_expr(const Point<dim> &p) const {
                    return 0.0;
                }
                double z_expr(const Point<dim> &p) const {
                    return 0.0;
                }
        };

        template<int dim>
        class NeumannCondition0 : public Function<dim> {
            public:
            NeumannCondition0(const unsigned int components) : Function<dim>(components) {}

            virtual
            void vector_value(const Point<dim> &p, Vector<double> &values) const override {
                values[0] = x_expr(p);
                if constexpr (dim >= 2) values[1] = y_expr(p);
                if constexpr (dim >= 3) values[2] = z_expr(p);
            }

            virtual
            double value(const Point<dim> &p, const unsigned int component = 0) const override {
                if (component == 0) return x_expr(p);
                if (component == 1) return y_expr(p);
                return z_expr(p);
            }

            double scalar_value() const {
                return 0.0;
            }

            protected:
            double x_expr(const Point<dim> &p) const {
                return 0.0; // -p_out [Pa]
            }
            double y_expr(const Point<dim> &p) const {
                return 0.0;
            }
            double z_expr(const Point<dim> &p) const {
                return 0.0;
            }
        };

        template<int dim>
        class NeumannCondition1 : public Function<dim> {
            public:
            NeumannCondition1(const unsigned int components) : Function<dim>(components) {}

            virtual
            void vector_value(const Point<dim> &p, Vector<double> &values) const override {
                values[0] = x_expr(p);
                if constexpr (dim >= 2) values[1] = y_expr(p);
                if constexpr (dim >= 3) values[2] = z_expr(p);
            }

            virtual
            double value(const Point<dim> &p, const unsigned int component = 0) const override {
                if (component == 0) return x_expr(p);
                if (component == 1) return y_expr(p);
                return z_expr(p);
            }

            double scalar_value() const {
                return 0.0;
            }

            protected:
            double x_expr(const Point<dim> &p) const {
                return 0.0;
            }
            double y_expr(const Point<dim> &p) const {
                return 0.0;
            }
            double z_expr(const Point<dim> &p) const {
                return 0.0;
            }
        };
    }

    namespace Preconditioners {

        constexpr unsigned int solver_iters = 1000;
        constexpr double tol = 1e-2;
        constexpr unsigned int vel_blk = 0;
        constexpr unsigned int press_blk = 1;

        /*
         * 1 0
         * 0 1
         */
        class Identity {
            public:
            void vmult(BlkVector &dst, const BlkVector &src) const {
                dst = src;
            }
        };

        /*
         * P =
         * A 0
         * 0 (1/nu) * M
         *
         * P^-1 =
         * A^-1 0
         * 0 nu*M^-1
         */
        class BlockDiagonal {
            protected:
            const SpMtx *vel_stiff;
            ILU prec_vel;
            const SpMtx *pressure_mass;
            ILU prec_pressure;

            public:
            void initialize(const SpMtx &vel_stiff_, const SpMtx &pressure_mass_) {
                vel_stiff = &vel_stiff_;
                pressure_mass = &pressure_mass_;
                prec_vel.initialize(vel_stiff_);
                prec_pressure.initialize(pressure_mass_);
            }

            void vmult(BlkVector &dst, const BlkVector &src) const {
                const auto vel_l2_norm = src.block(vel_blk).l2_norm();
                const auto press_l2_norm = src.block(press_blk).l2_norm();
                SolverControl control_vel(solver_iters,tol * vel_l2_norm);
                SolverControl control_press(solver_iters,tol * press_l2_norm);

                SolverCG<DVector> cg_vel(control_vel);
                SolverCG<DVector> cg_press(control_press);

                cg_vel.solve(
                    *vel_stiff,
                    dst.block(vel_blk),
                    src.block(vel_blk),
                    prec_vel
                );
                cg_press.solve(
                    *pressure_mass,
                    dst.block(press_blk),
                    src.block(press_blk),
                    prec_pressure
                );
            }


        };

        /*
         * P =
         * A 0
         * B (1/nu) * M
         *
         * P^-1 [x,y] = [A^-1x, nu*M^-1(y-BA^-1x)] = [u,p]
         */
        class BlockTriangular {
            protected:
            const SpMtx *vel_stiff;
            ILU prec_vel;
            const SpMtx *pressure_mass;
            ILU prec_pressure;
            const SpMtx *B;
            // temporary vector for results
            mutable DVector tmp;

            public:

            void initialize(
                const SpMtx &vel_stiff_,
                const SpMtx &pressure_mass_,
                const SpMtx &B_
                ) {
                vel_stiff = &vel_stiff_;
                pressure_mass = &pressure_mass_;
                B = &B_;
                prec_vel.initialize(vel_stiff_);
                prec_pressure.initialize(pressure_mass_);
            }

            void vmult(BlkVector &dst, const BlkVector &src) const {
                const auto vel_l2_norm = src.block(vel_blk).l2_norm();
                const auto press_l2_norm = src.block(press_blk).l2_norm();
                SolverControl control_vel(solver_iters,tol * vel_l2_norm);
                SolverControl control_press(solver_iters,tol * press_l2_norm);

                SolverCG<DVector> cg_vel(control_vel);
                SolverCG<DVector> cg_press(control_press);

                // u = A^-1 * x
                cg_vel.solve(
                    *vel_stiff,
                    dst.block(vel_blk),
                    src.block(vel_blk),
                    prec_vel
                );

                // tmp dimensions match p
                tmp.reinit(src.block(press_blk));
                // tmp = B * A^-1 x
                B->vmult(tmp,dst.block(vel_blk));
                // tmp = y - (B * A^-1 x)
                tmp.sadd(-1.0, src.block(press_blk));
                // p = M^-1 * (y - (B * A^-1 x))
                cg_press.solve(
                    *pressure_mass,
                    dst.block(press_blk),
                    tmp,
                    prec_pressure
                );
            }

        };

    }

}

#endif //NNPDE_PROJECT_NAVIERSTOKESHELPERS_HPP