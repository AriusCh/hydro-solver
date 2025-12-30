#ifndef HYDRO_SOVLER_PROBLEM_HPP_
#define HYDRO_SOVLER_PROBLEM_HPP_

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "eos.hpp"
#include "material_closure.hpp"
#include "plastic_properties.hpp"

enum class ProblemDimension { e1D, e2D };

enum class BoundaryType {
  eFree,
  eWall,
  eNoSlipWall,
  eTransmissive,
};

struct Problem {
  std::string name;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double tmin;
  double tmax;
  std::vector<double> tOut;
  double tMul;

  ProblemDimension dimension;

  BoundaryType leftBoundaryType;
  BoundaryType topBoundaryType;
  BoundaryType rightBoundaryType;
  BoundaryType bottomBoundaryType;

  std::size_t numberOfMaterials;

  std::function<double(const double x, const double y)> uInitializer;
  std::function<double(const double x, const double y)> vInitializer;
  std::function<std::vector<double>(const double x, const double y)>
      volFractionInitializer;
  std::function<std::vector<double>(const double x, const double y)>
      rhoInitializer;
  std::function<std::vector<double>(const double x, const double y)>
      pInitializer;

  std::vector<std::shared_ptr<EOS>> eoses;

  std::shared_ptr<MaterialClosure> matClosure;
  std::shared_ptr<PlasticProperties> plastic_params = nullptr;
};

struct DefaultProblems {
  static const Problem sodTest;
  static const Problem vacuumShock;
  static const Problem blastWave;
  static const Problem laserVolumeTarget;
  static const Problem laserVolumeTargetFreqOut;
  static const Problem plasticLaserVolumeTarget;
  static const Problem laserVolumeTargetPlasticSolid;
  static const Problem laserVolumeTargetPlasticSolidExtended;
  static const Problem laserVolumeTargetWithSeparation;
  static const Problem laserVolumeTargetWithSeparationFreqOut;
  static const Problem laserVolumeTargetWithCellSeparation;
  static const Problem laserVolumeTargetWithCellSeparationFreqOut;
  static const Problem twoMatSodTest;
  static const Problem waterTube;
  static const Problem liquidGas;
  static const Problem plasticImpact;
  static const Problem plasticImpactFreqOut;
  static const Problem plasticPistonLike;
  static const Problem plasticWilkins;
  static const Problem plasticRiemann;
  static const Problem plasticRiemannCorrect;
  static const Problem plasticNailing;
};

// factory functions
Problem createRiemannProblem(const std::string &name, const double xmin,
                             const double xmax, const double spl,
                             const double tmax, const double tMul,
                             const double uL, const double rhoL,
                             const double pL, const double uR,
                             const double rhoR, const double pR,
                             const double gamma);
Problem createMultiMaterialRiemannProblem(
    const std::string &name, const double xmin, const double xmax,
    const double spl, const double tmax, const double tMul, const double uL,
    const double rhoL, const double pL, const double gammaL, const double p0L,
    const double uR, const double rhoR, const double pR, const double gammaR,
    const double p0R);
Problem createVacuumShockProblem(const std::string &name, const double xmin,
                                 const double xmax, const double tmax,
                                 const double tMul, const double u,
                                 const double rho, const double p,
                                 const double gamma);
Problem createBlastWaveProblem(
    const std::string &name, const double xmin, const double xmax,
    const double ymin, const double ymax, const double spl, const double tmax,
    const double tMul, const double uIn, const double vIn, const double rhoIn,
    const double pIn, const double uOut, const double vOut, const double rhoOut,
    const double pOut, const double gamma);
Problem createLaserVolumeTargetProblem(const std::string &name, double xmin,
                                       double xmax, double ymin, double tmax,
                                       const std::vector<double> &tOut,
                                       double rhoM, double pCold, double pHeat,
                                       double RL, double dSkin);
Problem createLaserVolumeTargetPlasticSolidProblem(
    const std::string &name, double xmin, double xmax, double ymin, double tmax,
    const std::vector<double> &tOut, double rhoM, double pCold, double pHeat,
    double RL, double dSkin);
Problem createLaserVolumeTargetWithSeparationProblem(
    const std::string &name, double xmin, double xmax, double ymin, double tmax,
    const std::vector<double> &tOut, double rhoM, double pCold, double pHeat,
    double RL, double dSkin);
Problem createLaserVolumeTargetWithCellSeparationProblem(
    const std::string &name, double xmin, double xmax, double ymin, double tmax,
    const std::vector<double> &tOut, double rhoM, double pCold, double pHeat,
    double RL, double dSkin);

Problem createPlasticImpactProblem(const std::string &name, const double l,
                                   const double h, const double tmax,
                                   const std::vector<double> &tOut,
                                   const double tMul, const double u0);
Problem createPlasticPistonLike(const std::string &name, const double tmax,
                                const std::vector<double> &tOut,
                                const double tMul);
Problem createPlasticWilkins(const std::string &name, const double tmax,
                             const std::vector<double> &tOut,
                             const double tMul);
Problem createPlasticRiemann(const std::string &name, const double tmax,
                             const std::vector<double> &tOut, const double tMul,
                             const double bulk_modulus, const double rho0,
                             const double shear_modulus,
                             const double equivalent_stress, const double L);
Problem createPlasticRiemannCorrect(const std::string &name, const double tmax,
                                    const std::vector<double> &tOut,
                                    const double tMul);
Problem createPlasticNailingProblem(
    const std::string &name, const double xmax, const double ymax,
    const double x_split, const double y_split, const double tmax,
    const std::vector<double> &tOut, const double tMul, const double u_left,
    const double bulk_modulus, const double rho0, const double shear_modulus,
    const double equivalent_stress);

#endif  // HYDRO_SOVLER_PROBLEM_HPP_
