#include <Eigen/Core>

// #include "methods/plastic_elastic_fem_method.hpp"
#include "methods/lagrangian_fem_method.hpp"
#include "problem.hpp"
#include "simulation.hpp"

using MethodType = LagrangianFemMethod;

int main() {
  {
    const Problem &problem =
        createRiemannProblem("sod-test", 0.0, 1.0, 0.5, 0.3, 1.0, 0.0, 1.0, 1.0,
                             0.0, 0.125, 0.1, 5.0 / 3.0);
    std::unique_ptr<MethodType> mtd =
        std::make_unique<MethodType>(problem, 100, 1, 2);
    Simulation sim(std::move(mtd));
    sim.run();
  }
  {
    const Problem &problem =
        createRiemannProblem("sod-test", 0.0, 1.0, 0.5, 0.4, 1.0, 0.0, 1.0, 1.0,
                             0.0, 0.125, 0.1, 5.0 / 3.0);
    std::unique_ptr<MethodType> mtd =
        std::make_unique<MethodType>(problem, 100, 1, 2);
    Simulation sim(std::move(mtd));
    sim.run();
  }
  // {
  //   const Problem &problem =
  //   DefaultProblems::laserVolumeTargetPlasticSolidExtended;
  //   std::unique_ptr<MethodType> mtd =
  //       std::make_unique<MethodType>(problem, 1350, 1200, 1);
  //   Simulation sim(std::move(mtd));
  //   sim.run();
  // }
  // {
  //   const Problem &problem =
  //   DefaultProblems::laserVolumeTargetPlasticSolidExtended;
  //   std::unique_ptr<MethodType> mtd =
  //       std::make_unique<MethodType>(problem, 1350, 1200, 2);
  //   Simulation sim(std::move(mtd));
  //   sim.run();
  // }

  return 0;
}
