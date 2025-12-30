#ifndef HYDRO_SOLVER_MATERIAL_CLOSURES_PRESSURE_RELAX_CLOSURE_HPP_
#define HYDRO_SOLVER_MATERIAL_CLOSURES_PRESSURE_RELAX_CLOSURE_HPP_

#include "../material_closure.hpp"

class PressureRelaxClosure : public MaterialClosure {
 public:
  PressureRelaxClosure(const std::size_t kNumberOfMaterials_)
      : MaterialClosure(kNumberOfMaterials_) {}

  PressureRelaxClosure(const PressureRelaxClosure &) = delete;
  PressureRelaxClosure(PressureRelaxClosure &&) = delete;

  PressureRelaxClosure &operator=(const PressureRelaxClosure &) = delete;
  PressureRelaxClosure &operator=(PressureRelaxClosure &&) = delete;

  virtual ~PressureRelaxClosure() = default;

 public:
  virtual std::vector<double> calcVolFracRates(
      const std::vector<double> &volFracs, const std::vector<double> &rhos,
      const std::vector<double> &ps, const std::vector<double> &soundSpeeds,
      const std::size_t cell, const double h, const double dt);

  virtual double calcExchangePressure(const std::vector<double> &volFracs,
                                      const std::vector<double> &rhos,
                                      const std::vector<double> &ps,
                                      const std::vector<double> &soundSpeeds,
                                      const std::vector<double> &volFracRate,
                                      const double velocityScalarGrad,
                                      const double dt);
};

#endif  // HYDRO_SOLVER_MATERIAL_CLOSURES_PRESSURE_RELAX_CLOSURE_HPP_
