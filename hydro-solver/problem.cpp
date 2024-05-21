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
    std::vector<double> output(kNumberOfMatrials);
    if (x < spl) {
      output[0] = uL;
    } else {
      output[0] = uR;
    }
    return output;
  };
  auto vInitializer = []([[maybe_unused]] const double x,
                         [[maybe_unused]] const double y) {
    std::vector<double> output(kNumberOfMatrials, 0.0);
    return output;
  };
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
                 1,
                 uInitializer,
                 vInitializer,
                 volFractionInitializer,
                 rhoInitializer,
                 pInitializer,
                 eoses};
  return output;
}
