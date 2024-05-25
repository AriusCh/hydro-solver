#include <Eigen/Core>

#include "methods/lagrangian_fem_method.hpp"
#include "problem.hpp"
#include "simulation.hpp"

int main() {
  const Problem &problem = DefaultProblems::sodTest;
  std::unique_ptr<LagrangianFemMethod> mtd =
      std::make_unique<LagrangianFemMethod>(problem, 100, 1, 1);
  Simulation sim(std::move(mtd));
  sim.run();

  return 0;
}
