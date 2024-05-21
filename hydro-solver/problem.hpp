#ifndef HYDRO_SOVLER_PROBLEM_HPP_
#define HYDRO_SOVLER_PROBLEM_HPP_

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "eos.hpp"

enum class ProblemDimension { e1D, e2D };

enum class BoundaryType { eFree, eWall, eNoSlipWall };

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

  std::function<std::vector<double>(const double x, const double y)>
      uInitializer;
  std::function<std::vector<double>(const double x, const double y)>
      vInitializer;
  std::function<std::vector<double>(const double x, const double y)>
      volFractionInitializer;
  std::function<std::vector<double>(const double x, const double y)>
      rhoInitializer;
  std::function<std::vector<double>(const double x, const double y)>
      pInitializer;

  std::vector<std::shared_ptr<EOS>> eoses;
};

struct DefaultProblems {
  static const Problem sodTest;
};

// factory functions
Problem createRiemannProblem(const std::string &name, const double xmin,
                             const double xmax, const double spl,
                             const double tmax, const double tMul,
                             const double uL, const double rhoL,
                             const double pL, const double uR,
                             const double rhoR, const double pR,
                             const double gamma);

#endif  // HYDRO_SOVLER_PROBLEM_HPP_
