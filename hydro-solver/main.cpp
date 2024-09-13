#include <Eigen/Core>

#include "methods/lagrangian_fem_method.hpp"
#include "problem.hpp"
#include "simulation.hpp"

int main() {
  const Problem &problem = DefaultProblems::laserVolumeTargetWithSeparation;
  std::unique_ptr<LagrangianFemMethod> mtd =
      std::make_unique<LagrangianFemMethod>(problem, 90, 80, 2);
  Simulation sim(std::move(mtd));
  sim.run();

  return 0;
}
