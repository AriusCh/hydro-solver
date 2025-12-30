#include "problem.hpp"

#include <memory>

#include "eoses/eos_ideal_gas.hpp"
#include "eoses/eos_logarithmic.hpp"
#include "eoses/eos_mgal_precise6.hpp"
#include "eoses/eos_mgal_solid.hpp"
#include "eoses/eos_migruneisen_general.hpp"
#include "eoses/eos_simple_plastic.hpp"
#include "eoses/eos_stiffened_gas.hpp"
#include "material_closures/al_vac_cell_closure.hpp"
#include "material_closures/al_vac_closure.hpp"
#include "material_closures/pressure_relax_closure.hpp"
#include "plastic_properties.hpp"

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
                 BoundaryType::eTransmissive,
                 BoundaryType::eWall,
                 BoundaryType::eTransmissive,
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

Problem createMultiMaterialRiemannProblem(
    const std::string &name, const double xmin, const double xmax,
    const double spl, const double tmax, const double tMul, const double uL,
    const double rhoL, const double pL, const double gammaL, const double p0L,
    const double uR, const double rhoR, const double pR, const double gammaR,
    const double p0R) {
  constexpr std::size_t kNumberOfMatrials = 2;
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
  auto volFractionInitializer = [spl](const double x,
                                      [[maybe_unused]] const double y) {
    std::vector<double> output(kNumberOfMatrials);

    if (x < spl) {
      output[0] = 1.0;
      output[1] = 0.0;
    } else {
      output[0] = 0.0;
      output[1] = 1.0;
    }

    return output;
  };
  auto rhoInitializer = [rhoL, rhoR]([[maybe_unused]] const double x,
                                     [[maybe_unused]] const double y) {
    std::vector<double> output(kNumberOfMatrials);
    output[0] = rhoL;
    output[1] = rhoR;
    return output;
  };
  auto pInitializer = [pL, pR]([[maybe_unused]] const double x,
                               [[maybe_unused]] const double y) {
    std::vector<double> output(kNumberOfMatrials);
    output[0] = pL;
    output[1] = pR;
    return output;
  };
  std::vector<std::shared_ptr<EOS>> eoses{
      std::make_shared<EOSIdealGasWithBackgroundPressure>(gammaL, p0L),
      std::make_shared<EOSIdealGasWithBackgroundPressure>(gammaR, p0R)};
  assert(eoses.size() == kNumberOfMatrials);

  std::shared_ptr<MaterialClosure> matClosure =
      std::make_shared<PressureRelaxClosure>(kNumberOfMatrials);

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

  std::shared_ptr<PlasticProperties> plasticProperties =
      std::make_shared<PlasticProperties>();
  plasticProperties->shear_modulus = 80e9;
  plasticProperties->equivalent_stress = 12e9;

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
                 std::move(matClosure),
                 std::move(plasticProperties)};
  return output;
}

Problem createLaserVolumeTargetPlasticSolidProblem(
    const std::string &name, double xmin, double xmax, double ymin, double tmax,
    const std::vector<double> &tOut, double rhoM, double pCold, double pHeat,
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

  std::vector<std::shared_ptr<EOS>> eoses{std::make_shared<EOSMGAlSolid>()};

  std::shared_ptr<MaterialClosure> matClosure =
      std::make_shared<MaterialClosure>(kNumberOfMaterials);

  std::shared_ptr<PlasticProperties> plasticProperties =
      std::make_shared<PlasticProperties>();
  plasticProperties->shear_modulus = 60e9;
  plasticProperties->equivalent_stress = 12.6e9;

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
                 std::move(matClosure),
                 std::move(plasticProperties)};
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

Problem createLaserVolumeTargetWithCellSeparationProblem(
    const std::string &name, double xmin, double xmax, double ymin, double tmax,
    const std::vector<double> &tOut, double rhoM, double pCold, double pHeat,
    double RL, double dSkin) {
  Problem output = createLaserVolumeTargetProblem(
      name, xmin, xmax, ymin, tmax, tOut, rhoM, pCold, pHeat, RL, dSkin);
  output.matClosure = std::make_shared<AlVacCellClosure>();

  return output;
}

Problem createPlasticImpactProblem(const std::string &name, const double l,
                                   const double h, const double tmax,
                                   const std::vector<double> &tOut,
                                   const double tMul, const double u0) {
  constexpr std::size_t kNumberOfMatrials = 1;
  auto uInitializer = [u0](const double x, [[maybe_unused]] const double y) {
    if (x <= 0.0) {
      return 0.0;
    } else {
      return u0;
    }
  };
  auto vInitializer = []([[maybe_unused]] const double x,
                         [[maybe_unused]] const double y) { return 0.0; };

  auto volFractionInitializer = []([[maybe_unused]] const double x,
                                   [[maybe_unused]] const double y) {
    std::vector<double> output(kNumberOfMatrials, 1.0);
    return output;
  };

  auto rhoInitializer = []([[maybe_unused]] const double x,
                           [[maybe_unused]] const double y) {
    std::vector<double> output(kNumberOfMatrials);
    output[0] = 2785.0;
    return output;
  };
  auto pInitializer = []([[maybe_unused]] const double x,
                         [[maybe_unused]] const double y) {
    std::vector<double> output(kNumberOfMatrials);
    output[0] = 1e-6;
    return output;
  };

  std::vector<std::shared_ptr<EOS>> eoses{
      std::make_shared<EOSGruneisesenGeneral>(2785.0, 5328.0, 2.0, 1.338)};

  std::shared_ptr<MaterialClosure> matClosure =
      std::make_shared<MaterialClosure>(kNumberOfMatrials);

  std::shared_ptr<PlasticProperties> plasticProperties =
      std::make_shared<PlasticProperties>();
  plasticProperties->shear_modulus = 27.6e9;
  plasticProperties->equivalent_stress = 300e6;

  Problem output{name,
                 0.0,
                 l,
                 -0.5 * h,
                 0.5 * h,
                 0.0,
                 tmax,
                 tOut,
                 tMul,
                 ProblemDimension::e2D,
                 BoundaryType::eWall,
                 BoundaryType::eFree,
                 BoundaryType::eFree,
                 BoundaryType::eFree,
                 kNumberOfMatrials,
                 uInitializer,
                 vInitializer,
                 volFractionInitializer,
                 rhoInitializer,
                 pInitializer,
                 eoses,
                 std::move(matClosure),
                 std::move(plasticProperties)};
  return output;
}

Problem createPlasticPistonLike(const std::string &name, const double tmax,
                                const std::vector<double> &tOut,
                                const double tMul) {
  constexpr std::size_t kNumberOfMatrials = 1;
  auto uInitializer = [](const double x, [[maybe_unused]] const double y) {
    if (x <= 0.0) {
      return 20.0;
    } else {
      return 0.0;
    }
  };
  auto vInitializer = []([[maybe_unused]] const double x,
                         [[maybe_unused]] const double y) { return 0.0; };

  auto volFractionInitializer = []([[maybe_unused]] const double x,
                                   [[maybe_unused]] const double y) {
    std::vector<double> output(kNumberOfMatrials, 1.0);
    return output;
  };

  auto rhoInitializer = []([[maybe_unused]] const double x,
                           [[maybe_unused]] const double y) {
    std::vector<double> output(kNumberOfMatrials);
    output[0] = 8930.0;
    return output;
  };
  auto pInitializer = []([[maybe_unused]] const double x,
                         [[maybe_unused]] const double y) {
    std::vector<double> output(kNumberOfMatrials);
    output[0] = 1e5;
    return output;
  };

  std::vector<std::shared_ptr<EOS>> eoses{
      std::make_shared<EOSGruneisesenGeneral>(8930.0, 3940.0, 2.0, 1.49)};

  std::shared_ptr<MaterialClosure> matClosure =
      std::make_shared<MaterialClosure>(kNumberOfMatrials);

  std::shared_ptr<PlasticProperties> plasticProperties =
      std::make_shared<PlasticProperties>();
  plasticProperties->shear_modulus = 45e9;
  plasticProperties->equivalent_stress = 90e6;

  Problem output{name,
                 0.0,
                 1,
                 0.0,
                 0.1,
                 0.0,
                 tmax,
                 tOut,
                 tMul,
                 ProblemDimension::e2D,
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
                 std::move(matClosure),
                 std::move(plasticProperties)};
  return output;
}
Problem createPlasticWilkins(const std::string &name, const double tmax,
                             const std::vector<double> &tOut,
                             const double tMul) {
  constexpr std::size_t kNumberOfMaterials = 1;
  constexpr double x_split = 5e-3;

  auto uInitializer = [](const double x, [[maybe_unused]] const double y) {
    if (x < x_split) {
      return 800.0;
    } else {
      return 0.0;
    }
  };
  auto vInitializer = []([[maybe_unused]] const double x,
                         [[maybe_unused]] const double y) { return 0.0; };

  auto volFractionInitializer = []([[maybe_unused]] const double x,
                                   [[maybe_unused]] const double y) {
    std::vector<double> output(kNumberOfMaterials, 1.0);
    return output;
  };

  auto rhoInitializer = []([[maybe_unused]] const double x,
                           [[maybe_unused]] const double y) {
    std::vector<double> output(kNumberOfMaterials);
    output[0] = 2785.0;
    return output;
  };
  auto pInitializer = []([[maybe_unused]] const double x,
                         [[maybe_unused]] const double y) {
    std::vector<double> output(kNumberOfMaterials);
    output[0] = 1e-6;
    return output;
  };

  std::vector<std::shared_ptr<EOS>> eoses{
      std::make_shared<EOSGruneisesenGeneral>(2785.0, 5328.0, 2.0, 1.338)};

  std::shared_ptr<MaterialClosure> matClosure =
      std::make_shared<MaterialClosure>(kNumberOfMaterials);

  std::shared_ptr<PlasticProperties> plasticProperties =
      std::make_shared<PlasticProperties>();
  plasticProperties->shear_modulus = 27.6e9;
  plasticProperties->equivalent_stress = 300e6;

  Problem output{name,
                 0.0,
                 50e-3,
                 0.0,
                 5e-3,
                 0.0,
                 tmax,
                 tOut,
                 tMul,
                 ProblemDimension::e2D,
                 BoundaryType::eFree,
                 BoundaryType::eWall,
                 BoundaryType::eWall,
                 BoundaryType::eWall,
                 kNumberOfMaterials,
                 uInitializer,
                 vInitializer,
                 volFractionInitializer,
                 rhoInitializer,
                 pInitializer,
                 eoses,
                 std::move(matClosure),
                 std::move(plasticProperties)};
  return output;
}

Problem createPlasticRiemann(const std::string &name, const double tmax,
                             const std::vector<double> &tOut, const double tMul,
                             const double bulk_modulus, const double rho0,
                             const double shear_modulus,
                             const double equivalent_stress, const double L) {
  constexpr std::size_t kNumberOfMaterials = 1;
  const double x_split = 0.5 * L;

  constexpr double p_left = 4e9;
  const double rho_left = rho0 * (1.0 + p_left / bulk_modulus);

  auto uInitializer = []([[maybe_unused]] const double x,
                         [[maybe_unused]] const double y) { return 0.0; };
  auto vInitializer = []([[maybe_unused]] const double x,
                         [[maybe_unused]] const double y) { return 0.0; };

  auto volFractionInitializer = []([[maybe_unused]] const double x,
                                   [[maybe_unused]] const double y) {
    std::vector<double> output(kNumberOfMaterials, 1.0);
    return output;
  };

  auto rhoInitializer = [rho_left, rho0, x_split](
                            [[maybe_unused]] const double x,
                            [[maybe_unused]] const double y) {
    std::vector<double> output(kNumberOfMaterials);
    if (x < x_split) {
      output[0] = rho_left;
    } else {
      output[0] = rho0;
    }
    return output;
  };

  auto pInitializer = [x_split](const double x,
                                [[maybe_unused]] const double y) {
    std::vector<double> output(kNumberOfMaterials);
    if (x < x_split) {
      output[0] = p_left;
    } else {
      output[0] = 0.0e9;
    }
    return output;
  };

  std::vector<std::shared_ptr<EOS>> eoses{
      std::make_shared<EOSSimplePlastic>(bulk_modulus, rho0)};

  std::shared_ptr<MaterialClosure> matClosure =
      std::make_shared<MaterialClosure>(kNumberOfMaterials);

  std::shared_ptr<PlasticProperties> plasticProperties =
      std::make_shared<PlasticProperties>();
  plasticProperties->shear_modulus = shear_modulus;
  plasticProperties->equivalent_stress = equivalent_stress;

  Problem output{name,
                 0.0,
                 L,
                 0.0,
                 0.1 * L,
                 0.0,
                 tmax,
                 tOut,
                 tMul,
                 ProblemDimension::e1D,
                 BoundaryType::eWall,
                 BoundaryType::eWall,
                 BoundaryType::eWall,
                 BoundaryType::eWall,
                 kNumberOfMaterials,
                 uInitializer,
                 vInitializer,
                 volFractionInitializer,
                 rhoInitializer,
                 pInitializer,
                 eoses,
                 std::move(matClosure),
                 std::move(plasticProperties)};
  return output;
}

Problem createPlasticRiemannCorrect(const std::string &name, const double tmax,
                                    const std::vector<double> &tOut,
                                    const double tMul) {
  constexpr std::size_t kNumberOfMaterials = 1;
  constexpr double x_split = 0.0;

  constexpr double Y0 = 0.3e9;
  constexpr double rho0 = 2688.9;
  constexpr double p_left = 3.8e9;
  constexpr double B = 73e9;
  constexpr double G = 23e9;
  const double rho_left = rho0 * std::exp(p_left / B);

  auto uInitializer = []([[maybe_unused]] const double x,
                         [[maybe_unused]] const double y) { return 0.0; };
  auto vInitializer = []([[maybe_unused]] const double x,
                         [[maybe_unused]] const double y) { return 0.0; };

  auto volFractionInitializer = []([[maybe_unused]] const double x,
                                   [[maybe_unused]] const double y) {
    std::vector<double> output(kNumberOfMaterials, 1.0);
    return output;
  };

  auto rhoInitializer = [rho_left]([[maybe_unused]] const double x,
                                   [[maybe_unused]] const double y) {
    std::vector<double> output(kNumberOfMaterials);
    if (x < x_split) {
      output[0] = rho_left;
    } else {
      output[0] = rho0;
    }
    return output;
  };

  auto pInitializer = [](const double x, [[maybe_unused]] const double y) {
    std::vector<double> output(kNumberOfMaterials);
    if (x < x_split) {
      output[0] = p_left;
    } else {
      output[0] = 0.0e9;
    }
    return output;
  };

  std::vector<std::shared_ptr<EOS>> eoses{std::make_shared<EOSLogarithmic>()};

  std::shared_ptr<MaterialClosure> matClosure =
      std::make_shared<MaterialClosure>(kNumberOfMaterials);

  std::shared_ptr<PlasticProperties> plasticProperties =
      std::make_shared<PlasticProperties>();
  plasticProperties->shear_modulus = G;
  plasticProperties->equivalent_stress = Y0;

  Problem output{name,
                 -0.02,
                 0.02,
                 0.0,
                 1.0,
                 0.0,
                 tmax,
                 tOut,
                 tMul,
                 ProblemDimension::e1D,
                 BoundaryType::eWall,
                 BoundaryType::eWall,
                 BoundaryType::eWall,
                 BoundaryType::eWall,
                 kNumberOfMaterials,
                 uInitializer,
                 vInitializer,
                 volFractionInitializer,
                 rhoInitializer,
                 pInitializer,
                 eoses,
                 std::move(matClosure),
                 std::move(plasticProperties)};
  return output;
}

Problem createPlasticNailingProblem(
    const std::string &name, const double xmax, const double ymax,
    const double x_split, const double y_split, const double tmax,
    const std::vector<double> &tOut, const double tMul, const double u_left,
    const double bulk_modulus, const double rho0, const double shear_modulus,
    const double equivalent_stress) {
  constexpr std::size_t kNumberOfMaterials = 1;

  auto uInitializer = [u_left, x_split, y_split](const double x,
                                                 const double y) {
    if (x <= x_split && std::abs(y - 2.0) <= y_split) {
      return u_left;
    }
    return 0.0;
  };
  auto vInitializer = []([[maybe_unused]] const double x,
                         [[maybe_unused]] const double y) { return 0.0; };
  auto volFractionInitializer = []([[maybe_unused]] const double x,
                                   [[maybe_unused]] const double y) {
    std::vector<double> output(kNumberOfMaterials, 1.0);
    return output;
  };
  auto rhoInitializer = []([[maybe_unused]] const double x,
                           [[maybe_unused]] const double y) {
    std::vector<double> output(kNumberOfMaterials);
    output[0] = 2700;
    return output;
  };

  auto pInitializer = []([[maybe_unused]] const double x,
                         [[maybe_unused]] const double y) {
    std::vector<double> output(kNumberOfMaterials);
    output[0] = 0.1 * 1e6;
    return output;
  };

  std::vector<std::shared_ptr<EOS>> eoses{
      std::make_shared<EOSStiffenedGas>(rho0, 5380, 2.67, bulk_modulus)};

  std::shared_ptr<MaterialClosure> matClosure =
      std::make_shared<MaterialClosure>(kNumberOfMaterials);

  std::shared_ptr<PlasticProperties> plasticProperties =
      std::make_shared<PlasticProperties>();
  plasticProperties->shear_modulus = shear_modulus;
  plasticProperties->equivalent_stress = equivalent_stress;

  Problem output{name,
                 -1.0,
                 xmax,
                 0.0,
                 2.0 * ymax,
                 0.0,
                 tmax,
                 tOut,
                 tMul,
                 ProblemDimension::e1D,
                 BoundaryType::eWall,
                 BoundaryType::eWall,
                 BoundaryType::eWall,
                 BoundaryType::eWall,
                 kNumberOfMaterials,
                 uInitializer,
                 vInitializer,
                 volFractionInitializer,
                 rhoInitializer,
                 pInitializer,
                 eoses,
                 std::move(matClosure),
                 std::move(plasticProperties)};
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
const Problem DefaultProblems::plasticLaserVolumeTarget =
    createLaserVolumeTargetProblem(
        "laser-al-plastic", 0.0, 900e-9, -800e-9, 115.2e-12,
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
const Problem DefaultProblems::laserVolumeTargetPlasticSolid =
    createLaserVolumeTargetPlasticSolidProblem(
        "laser-al-plastic-solid", 0.0, 900e-9, -800e-9, 115.2e-12,
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
        2730, 0.0, 47.2e9, 200e-9, 80e-9);
const Problem DefaultProblems::laserVolumeTargetPlasticSolidExtended =
    createLaserVolumeTargetPlasticSolidProblem(
        "laser-al-plastic-solid", 0.0, 2700e-9, -2400e-9, 384.0e-12,
        []() {
          constexpr double tmax = 384.0e-12;
          constexpr std::size_t iMax = 400;
          constexpr double dt = tmax / iMax;
          std::vector<double> output;

          for (std::size_t i = 1; i < iMax; i++) {
            output.push_back(i * dt);
          }

          return output;
        }(),
        2730, 0.0, 47.2e9, 200e-9, 80e-9);
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
const Problem DefaultProblems::laserVolumeTargetWithCellSeparation =
    createLaserVolumeTargetWithCellSeparationProblem(
        "laser-al-cell-sep", 0.0, 900e-9, -800e-9, 115.2e-12,
        std::vector<double>{0.1 * 9.6e-12, 1.0 * 9.6e-12, 2.0 * 9.6e-12,
                            3.0 * 9.6e-12, 4.0 * 9.6e-12, 5.0 * 9.6e-12,
                            6.0 * 9.6e-12, 7.0 * 9.6e-12, 8.0 * 9.6e-12,
                            9.0 * 9.6e-12, 10.0 * 9.6e-12, 11.0 * 9.6e-12},
        2413.0, 0.0, 35.6e9, 200e-9, 80e-9);
const Problem DefaultProblems::laserVolumeTargetWithCellSeparationFreqOut =
    createLaserVolumeTargetWithCellSeparationProblem(
        "laser-al-cell-sep", 0.0, 900e-9, -800e-9, 115.2e-12,
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
const Problem DefaultProblems::twoMatSodTest =
    createMultiMaterialRiemannProblem("two-mat-sod-test", 0.0, 1.0, 0.505, 0.2,
                                      1.0, 0.0, 1.0, 2.0, 2.0, 0.0, 0.0, 0.125,
                                      0.1, 1.4, 0.0);
const Problem DefaultProblems::waterTube = createMultiMaterialRiemannProblem(
    "water-tube", 0.0, 1.0, 0.701, 2.2e-4, 1e3, 0.0, 1e3, 1e9, 4.4, 6e8, 0.0,
    50.0, 1e5, 1.4, 0.0);
const Problem DefaultProblems::liquidGas = createMultiMaterialRiemannProblem(
    "liquid-gas", 0.0, 1.0, 0.8025, 291e-6, 1e6, 0.0, 1e3, 2e8, 4.4, 6e8, 0.0,
    50, 1e5, 1.4, 0.0);
const Problem DefaultProblems::plasticImpact =
    createPlasticImpactProblem("plastic-impact", 5.0, 1.0, 5e-3, {}, 1e3, -150);
const Problem DefaultProblems::plasticImpactFreqOut =
    createPlasticImpactProblem(
        "plastic-impact", 5.0, 1.0, 5e-3,
        []() {
          constexpr double tmax = 5e-3;
          constexpr std::size_t iMax = 50;
          constexpr double dt = tmax / iMax;
          std::vector<double> output;

          for (std::size_t i = 1; i < iMax; i++) {
            output.push_back(i * dt);
          }

          return output;
        }(),
        1e3, -150);
const Problem DefaultProblems::plasticPistonLike = createPlasticPistonLike(
    "plastic-piston-like", 150e-6,
    []() {
      constexpr double tmax = 150e-6;
      constexpr std::size_t iMax = 150;
      constexpr double dt = tmax / iMax;
      std::vector<double> output;

      for (std::size_t i = 1; i < iMax; i++) {
        output.push_back(i * dt);
      }

      return output;
    }(),
    1e6);
const Problem DefaultProblems::plasticWilkins = createPlasticWilkins(
    "plastic-wilkins", 5e-6,
    []() {
      constexpr double tmax = 5e-6;
      constexpr std::size_t iMax = 10;
      constexpr double dt = tmax / iMax;
      std::vector<double> output;

      for (std::size_t i = 1; i < iMax; i++) {
        output.push_back(i * dt);
      }

      return output;
    }(),
    1e6);
const Problem DefaultProblems::plasticRiemann = createPlasticRiemann(
    "plastic-riemann", 5e-6, {}, 1e6, 73e9, 2700, 23e9, 0.3e9, 0.1);
const Problem DefaultProblems::plasticRiemannCorrect =
    createPlasticRiemannCorrect("plastic-riemann-correct", 2e-6, {}, 1e6);
const Problem DefaultProblems::plasticNailing = createPlasticNailingProblem(
    "plastic-nailing", 4.0, 2.0, 2.0, 1.0, 0.15e-3,
    []() {
      constexpr double tmax = 0.15e-3;
      constexpr std::size_t iMax = 15;
      constexpr double dt = tmax / iMax;
      std::vector<double> output;

      for (std::size_t i = 1; i < iMax; i++) {
        output.push_back(i * dt);
      }

      return output;
    }(),
    10e3, 100, 74e9, 2710, 26.5e9, 0.3e9);
