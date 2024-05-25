#include "problem.hpp"

#include "eoses/eos_ideal_gas.hpp"

Problem createRiemannProblem(const std::string &name, const double xmin,
                             const double xmax, const double spl,
                             const double tmax, const double tMul,
                             const double uL, const double rhoL,
                             const double pL, const double uR,
                             const double rhoR, const double pR,
                             const double gamma) {
  constexpr std::size_t kNumberOfMatrials = 1;
  auto uInitializer = [spl, uL, uR](const double x,
                                    [[maybe_unused]] const double y) {
    double output;
    if (x < spl) {
      output = uL;
    } else {
      output = uR;
    }
    return output;
  };
  auto vInitializer = []([[maybe_unused]] const double x,
                         [[maybe_unused]] const double y) { return 0.0; };
  auto volFractionInitializer = []([[maybe_unused]] const double x,
                                   [[maybe_unused]] const double y) {
    std::vector<double> output(kNumberOfMatrials, 1.0);
    return output;
  };
  auto rhoInitializer = [spl, rhoL, rhoR](const double x,
                                          [[maybe_unused]] const double y) {
    std::vector<double> output(kNumberOfMatrials);
    if (x < spl) {
      output[0] = rhoL;

    } else {
      output[0] = rhoR;
    }
    return output;
  };
  auto pInitializer = [spl, pL, pR](const double x,
                                    [[maybe_unused]] const double y) {
    std::vector<double> output(kNumberOfMatrials);
    if (x < spl) {
      output[0] = pL;
    } else {
      output[0] = pR;
    }
    return output;
  };
  std::vector<std::shared_ptr<EOS>> eoses{std::make_shared<EOSIdealGas>(gamma)};

  std::shared_ptr<MaterialClosure> matClosure =
      std::make_shared<MaterialClosure>(kNumberOfMatrials);

  Problem output{name,
                 xmin,
                 xmax,
                 0.0,
                 1.0,
                 0.0,
                 tmax,
                 {},
                 tMul,
                 ProblemDimension::e1D,
                 BoundaryType::eWall,
                 BoundaryType::eWall,
                 BoundaryType::eWall,
                 BoundaryType::eWall,
                 kNumberOfMatrials,
                 uInitializer,
                 vInitializer,
                 volFractionInitializer,
                 rhoInitializer,
                 pInitializer,
                 eoses,
                 std::move(matClosure)};
  return output;
}

const Problem DefaultProblems::sodTest =
    createRiemannProblem("sod-test", 0.0, 1.0, 0.5, 0.2, 1.0, 0.0, 1.0, 1.0,
                         0.0, 0.125, 0.1, 5.0 / 3.0);
