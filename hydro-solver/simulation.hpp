#ifndef HYDRO_SOLVER_SRC_SIMULATION_HPP_
#define HYDRO_SOLVER_SRC_SIMULATION_HPP_

#include <memory>
#include <queue>

#include "method.hpp"

class Simulation {
 public:
  Simulation(std::unique_ptr<Method> method_) : method(std::move(method_)) {}

  Simulation(const Simulation& rhs) = delete;
  Simulation(Simulation&& rhs) = default;

  Simulation& operator=(const Simulation& rhs) = delete;
  Simulation& operator=(Simulation&& rhs) = default;

  ~Simulation() = default;

 public:
  void run();

 private:
  void logSimulationStart() const;
  void logSimulationEnd(const double simulationTime) const;
  void logIteration(std::size_t iterationNumber, double t, double dt,
                    double calcTime, double remTime) const;

 private:
  std::unique_ptr<Method> method;

  static constexpr std::size_t stepTimesMaxSize = 50;
  std::queue<double> stepTimes;
  double stepTimesSum = 0.0;
};

#endif  // HYDRO_SOLVER_SRC_SIMULATION_HPP_
