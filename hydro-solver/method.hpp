#ifndef HYDRO_SOLVER_METHOD_HPP_
#define HYDRO_SOLVER_METHOD_HPP_

#include <filesystem>
#include <string>

#include "problem.hpp"

class Method {
 public:
  Method(const std::string& name_, const Problem& problem_)
      : name(name_),
        problem(problem_),
        t(problem_.tmin),
        dt(problem_.tmax - problem_.tmin),
        outputDirPath(std::filesystem::path("output") / problem.name / name) {
    std::filesystem::create_directories(outputDirPath);
  }

  Method(const Method& rhs) = default;
  Method(Method&& rhs) = default;

  Method& operator=(const Method& rhs) = default;
  Method& operator=(Method&& rhs) = default;

  virtual ~Method() = default;

 public:
  virtual void dumpSolverInfo() const = 0;
  virtual void dumpData() const = 0;

  virtual void calcStep() = 0;

 public:
  std::string name;
  Problem problem;

  double t;
  double dt;

  std::filesystem::path outputDirPath;
};

#endif  // HYDRO_SOLVER_METHOD_HPP_
