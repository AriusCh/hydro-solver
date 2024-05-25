#include "simulation.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <deque>
#include <iostream>
#include <string>

void Simulation::run() {
  std::deque<double> tOut;
  {
    std::vector<double> tOutSorted = method->problem.tOut;
    std::sort(tOutSorted.begin(), tOutSorted.end());
    tOut = std::deque(tOutSorted.cbegin(), tOutSorted.cend());
  }

  std::size_t iterationNum = 0;
  logSimulationStart();

  method->dumpSolverInfo();
  method->dumpData();

  auto simStartTimePoint = std::chrono::high_resolution_clock::now();

  while (method->t < method->problem.tmax) {
    const double tmpDt = method->dt;
    if (!tOut.empty()) {
      const double tNext = tOut.front();
      if (method->t + method->dt > tNext) {
        method->dt = tNext - method->t;
      }
    }
    if (method->t + method->dt > method->problem.tmax) {
      method->dt = method->problem.tmax - method->t;
    }

    auto stepStartTimePoint = std::chrono::high_resolution_clock::now();

    method->calcStep();

    auto stepEndTimePoint = std::chrono::high_resolution_clock::now();

    const double calcStepTime =
        std::chrono::duration<double>(stepEndTimePoint - stepStartTimePoint)
            .count();

    stepTimes.push(calcStepTime);
    stepTimesSum += calcStepTime;
    if (stepTimes.size() > stepTimesMaxSize) {
      stepTimesSum -= stepTimes.front();
      stepTimes.pop();
    }
    const double remainingTime = (method->problem.tmax - method->t) /
                                 method->dt * stepTimesSum / stepTimes.size();

    logIteration(++iterationNum, method->t, method->dt, calcStepTime,
                 remainingTime);

    if (!tOut.empty()) {
      const double tNext = tOut.front();
      if (method->t >= tNext) {
        method->dumpData();
        tOut.pop_front();
        method->dt = tmpDt;
      }
    }
  }

  auto simEndTimePoint = std::chrono::high_resolution_clock::now();

  const double simulationTime =
      std::chrono::duration<double>(simEndTimePoint - simStartTimePoint)
          .count();
  logSimulationEnd(simulationTime);

  method->dumpData();
}

void Simulation::logSimulationStart() const {
  const std::string message = "PROBLEM: " + method->problem.name +
                              " METHOD: " + method->name +
                              " STARTING SIMULATION";
  std::cout << message << std::endl;
}

void Simulation::logSimulationEnd(const double simulationTime) const {
  std::cout << "PROBLEM: " << method->problem.name
            << " METHOD: " << method->name
            << " SIMULATION COMPLETE. TIME: " << simulationTime << std::endl;
}

void Simulation::logIteration(std::size_t iterationNumber, double t, double dt,
                              double calcTime, double remTime) const {
  std::cout << "PROBLEM: " << method->problem.name
            << " METHOD: " << method->name << " ITERATION: " << std::setw(7)
            << iterationNumber << " t: " << t << " dt: " << dt
            << " STEP TIME: " << calcTime << " remTime: " << std::floor(remTime)
            << std::endl;
  ;
}
