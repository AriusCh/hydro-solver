#include "problem.hpp"

#include "eoses/eos_ideal_gas.hpp"
#include "eoses/eos_mgal_precise6.hpp"
#include "material_closures/al_vac_closure.hpp"

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

Problem createVacuumShockProblem(const std::string &name, const double xmin,
                                 const double xmax, const double tmax,
                                 const double tMul, const double u,
                                 const double rho, const double p,
                                 const double gamma) {
  constexpr std::size_t kNumberOfMatrials = 1;

  auto uInitializer = [u]([[maybe_unused]] const double x,
                          [[maybe_unused]] const double y) { return u; };

  auto vInitializer = []([[maybe_unused]] const double x,
                         [[maybe_unused]] const double y) { return 0.0; };

  auto volFractionInitializer = []([[maybe_unused]] const double x,
                                   [[maybe_unused]] const double y) {
    std::vector<double> output(kNumberOfMatrials, 1.0);
    return output;
  };

  auto rhoInitializer = [rho]([[maybe_unused]] const double x,
                              [[maybe_unused]] const double y) {
    std::vector<double> output(kNumberOfMatrials);
    output[0] = rho;
    return output;
  };
  auto pInitializer = [p]([[maybe_unused]] const double x,
                          [[maybe_unused]] const double y) {
    std::vector<double> output(kNumberOfMatrials);
    output[0] = p;
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
                 BoundaryType::eFree,
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
                 std::move(matClosure)

  };

  return output;
}

Problem createBlastWaveProblem(
    const std::string &name, const double xmin, const double xmax,
    const double ymin, const double ymax, const double spl, const double tmax,
    const double tMul, const double uIn, const double vIn, const double rhoIn,
    const double pIn, const double uOut, const double vOut, const double rhoOut,
    const double pOut, const double gamma) {
  constexpr std::size_t kNumberOfMatrials = 1;

  auto uInitializer = [spl, uIn, uOut](const double x, const double y) {
    if (x * x + y * y <= spl * spl) {
      return uIn;
    }
    return uOut;
  };

  auto vInitializer = [spl, vIn, vOut](const double x, const double y) {
    if (x * x + y * y <= spl * spl) {
      return vIn;
    }
    return vOut;
  };

  auto volFracInitializer = []([[maybe_unused]] const double x,
                               [[maybe_unused]] const double y) {
    std::vector<double> output(kNumberOfMatrials, 0.0);
    output[0] = 1.0;
    return output;
  };

  auto rhoInitializer = [spl, rhoIn, rhoOut](const double x, const double y) {
    std::vector<double> output(kNumberOfMatrials, 0.0);
    if (x * x + y * y <= spl * spl) {
      output[0] = rhoIn;
    } else {
      output[0] = rhoOut;
    }
    return output;
  };

  auto pInitializer = [spl, pIn, pOut](const double x, const double y) {
    std::vector<double> output(kNumberOfMatrials, 0.0);
    if (x * x + y * y <= spl * spl) {
      output[0] = pIn;
    } else {
      output[0] = pOut;
    }
    return output;
  };

  std::vector<std::shared_ptr<EOS>> eoses{std::make_shared<EOSIdealGas>(gamma)};

  std::shared_ptr<MaterialClosure> matClosure =
      std::make_shared<MaterialClosure>(kNumberOfMatrials);

  Problem output{name,
                 xmin,
                 xmax,
                 ymin,
                 ymax,
                 0.0,
                 tmax,
                 {},
                 tMul,
                 ProblemDimension::e2D,
                 BoundaryType::eWall,
                 BoundaryType::eWall,
                 BoundaryType::eWall,
                 BoundaryType::eWall,
                 kNumberOfMatrials,
                 uInitializer,
                 vInitializer,
                 volFracInitializer,
                 rhoInitializer,
                 pInitializer,
                 eoses,
                 std::move(matClosure)};
  return output;
}

Problem createLaserVolumeTargetProblem(const std::string &name, double xmin,
                                       double xmax, double ymin, double tmax,
                                       const std::vector<double> &tOut,
                                       double rhoM, double pCold, double pHeat,
                                       double RL, double dSkin) {
  constexpr std::size_t kNumberOfMaterials = 1;

  const double ymax = 0.0;
  const double tmin = 0.0;
  const double tMul = 1e12;
  BoundaryType leftBoundaryType = BoundaryType::eWall;
  BoundaryType topBoundaryType = BoundaryType::eFree;
  BoundaryType rightBoundaryType = BoundaryType::eWall;
  BoundaryType bottomBoundaryType = BoundaryType::eWall;
  ProblemDimension dimension = ProblemDimension::e2D;

  auto uInitializer = []([[maybe_unused]] double x, [[maybe_unused]] double y) {
    return 0.0;
  };
  auto vInitializer = []([[maybe_unused]] double x, [[maybe_unused]] double y) {
    return 0.0;
  };

  auto volFracInitializer = []([[maybe_unused]] double x,
                               [[maybe_unused]] double y) {
    std::vector<double> output(kNumberOfMaterials, 0.0);
    output[0] = 1.0;
    return output;
  };

  auto rhoInitializer = [rhoM]([[maybe_unused]] double x,
                               [[maybe_unused]] double y) {
    std::vector<double> output(kNumberOfMaterials, 0.0);
    output[0] = rhoM;
    return output;
  };

  auto pInitializer = [pCold, pHeat, RL, dSkin](double x, double y) {
    std::vector<double> output(kNumberOfMaterials, 0.0);
    if (std::abs(x) <= RL && y > -dSkin) {
      output[0] = pHeat;
    } else {
      output[0] = pCold;
    }
    return output;
  };

  std::vector<std::shared_ptr<EOS>> eoses{std::make_shared<EOSMGAlPrecise6>()};

  std::shared_ptr<MaterialClosure> matClosure =
      std::make_shared<MaterialClosure>(kNumberOfMaterials);

  Problem output{name,
                 xmin,
                 xmax,
                 ymin,
                 ymax,
                 tmin,
                 tmax,
                 tOut,
                 tMul,
                 dimension,
                 leftBoundaryType,
                 topBoundaryType,
                 rightBoundaryType,
                 bottomBoundaryType,
                 kNumberOfMaterials,
                 uInitializer,
                 vInitializer,
                 volFracInitializer,
                 rhoInitializer,
                 pInitializer,
                 eoses,
                 std::move(matClosure)};
  return output;
}

Problem createLaserVolumeTargetWithSeparationProblem(
    const std::string &name, double xmin, double xmax, double ymin, double tmax,
    const std::vector<double> &tOut, double rhoM, double pCold, double pHeat,
    double RL, double dSkin) {
  Problem output = createLaserVolumeTargetProblem(
      name, xmin, xmax, ymin, tmax, tOut, rhoM, pCold, pHeat, RL, dSkin);
  output.matClosure = std::make_shared<AlVacClosure>();

  return output;
}

const Problem DefaultProblems::sodTest =
    createRiemannProblem("sod-test", 0.0, 1.0, 0.5, 0.2, 1.0, 0.0, 1.0, 1.0,
                         0.0, 0.125, 0.1, 5.0 / 3.0);
const Problem DefaultProblems::vacuumShock = createVacuumShockProblem(
    "vacuum-shock", 0.0, 1.0, 0.2, 1.0, 0.0, 1.0, 1.0, 1.4);
const Problem DefaultProblems::blastWave =
    createBlastWaveProblem("blast-wave", 0.0, 1.0, 0.0, 1.0, 0.4, 0.25, 1.0,
                           0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.125, 0.1, 1.4);
const Problem DefaultProblems::laserVolumeTarget =
    createLaserVolumeTargetProblem(
        "laser-al", 0.0, 900e-9, -800e-9, 115.2e-12,
        std::vector<double>{0.1 * 9.6e-12, 1.0 * 9.6e-12, 2.0 * 9.6e-12,
                            3.0 * 9.6e-12, 4.0 * 9.6e-12, 5.0 * 9.6e-12,
                            6.0 * 9.6e-12, 7.0 * 9.6e-12, 8.0 * 9.6e-12,
                            9.0 * 9.6e-12, 10.0 * 9.6e-12, 11.0 * 9.6e-12},
        2413.0, 0.0, 35.6e9, 200e-9, 80e-9);
const Problem DefaultProblems::laserVolumeTargetFreqOut =
    createLaserVolumeTargetProblem(
        "laser-al", 0.0, 900e-9, -800e-9, 115.2e-12,
        []() {
          constexpr double tmax = 115.2e-12;
          constexpr std::size_t iMax = 120;
          constexpr double dt = tmax / iMax;
          std::vector<double> output;

          for (std::size_t i = 1; i < iMax; i++) {
            output.push_back(i * dt);
          }

          return output;
        }(),
        2413.0, 0.0, 35.6e9, 200e-9, 80e-9);
const Problem DefaultProblems::laserVolumeTargetWithSeparation =
    createLaserVolumeTargetWithSeparationProblem(
        "laser-al-sep", 0.0, 900e-9, -800e-9, 115.2e-12,
        std::vector<double>{0.1 * 9.6e-12, 1.0 * 9.6e-12, 2.0 * 9.6e-12,
                            3.0 * 9.6e-12, 4.0 * 9.6e-12, 5.0 * 9.6e-12,
                            6.0 * 9.6e-12, 7.0 * 9.6e-12, 8.0 * 9.6e-12,
                            9.0 * 9.6e-12, 10.0 * 9.6e-12, 11.0 * 9.6e-12},
        2413.0, 0.0, 35.6e9, 200e-9, 80e-9);
const Problem DefaultProblems::laserVolumeTargetWithSeparationFreqOut =
    createLaserVolumeTargetWithSeparationProblem(
        "laser-al-sep", 0.0, 900e-9, -800e-9, 115.2e-12,
        []() {
          constexpr double tmax = 115.2e-12;
          constexpr std::size_t iMax = 120;
          constexpr double dt = tmax / iMax;
          std::vector<double> output;

          for (std::size_t i = 1; i < iMax; i++) {
            output.push_back(i * dt);
          }

          return output;
        }(),
        2413.0, 0.0, 35.6e9, 200e-9, 80e-9);
