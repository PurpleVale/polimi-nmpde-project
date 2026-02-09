#include "Heat.hpp"

// Main function.
int
main(int argc, char *argv[])
{
  constexpr unsigned int dim = Heat::dim;

  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  const auto mu = [](const Point<dim> & /*p*/) { return 0.1; };
  const auto f  = [](const Point<dim>  &/*p*/, const double  &/*t*/) {
    return 0.0;
  };

  Heat problem(/*mesh_filename = */ "./text/meshes/mesh-cube-10.msh",
               /* degree = */ 1,
               /* T = */ 1.0,
               /* theta = */ 1.0,
               /* delta_t = */ 0.05,
               mu,
               f);

  problem.run();

  return 0;
}