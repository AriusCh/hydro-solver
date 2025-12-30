#ifndef HYDRO_SOLVER_MATERIAL_CLOSURE_HPP_
#define HYDRO_SOLVER_MATERIAL_CLOSURE_HPP_

#include <cassert>
#include <cstddef>
#include <vector>

class MaterialClosure {
 public:
  MaterialClosure(const std::size_t kNumberOfMaterials_)
      : kNumberOfMaterials(kNumberOfMaterials_) {}

  MaterialClosure(const MaterialClosure &) = delete;
  MaterialClosure(MaterialClosure &&) = delete;

  MaterialClosure &operator=(const MaterialClosure &) = delete;
  MaterialClosure &operator=(MaterialClosure &&) = delete;

  virtual ~MaterialClosure() = default;

 public:
  virtual std::vector<double> calcVolFracRates(
      const std::vector<double> &volFracs, const std::vector<double> &rhos,
      const std::vector<double> &ps, const std::vector<double> &soundSpeeds,
      const std::size_t cell, const double h, const double dt);
  virtual double calcExchangePressure(const std::vector<double> &volFracs,
                                      const std::vector<double> &rhos,
                                      const std::vector<double> &ps,
                                      const std::vector<double> &soundSpeeds,
                                      const std::vector<double> &volFracRates,
                                      const double velocityScalarGrad,
                                      const double dt);

 protected:
  const std::size_t kNumberOfMaterials;
};

inline std::vector<double> MaterialClosure::calcVolFracRates(
    [[maybe_unused]] const std::vector<double> &volFracs,
    [[maybe_unused]] const std::vector<double> &rhos,
    [[maybe_unused]] const std::vector<double> &ps,
    [[maybe_unused]] const std::vector<double> &soundSpeeds,
    [[maybe_unused]] std::size_t cell, [[maybe_unused]] const double h,
    [[maybe_unused]] const double dt) {
  assert(volFracs.size() == kNumberOfMaterials);
  assert(rhos.size() == kNumberOfMaterials);
  assert(ps.size() == kNumberOfMaterials);
  assert(soundSpeeds.size() == kNumberOfMaterials);
  return std::vector<double>(kNumberOfMaterials, 0.0);
}

inline double MaterialClosure::calcExchangePressure(
    [[maybe_unused]] const std::vector<double> &volFracs,
    [[maybe_unused]] const std::vector<double> &rhos,
    [[maybe_unused]] const std::vector<double> &ps,
    [[maybe_unused]] const std::vector<double> &soundSpeeds,
    [[maybe_unused]] const std::vector<double> &volFracRates,
    [[maybe_unused]] const double velocityScalarGrad,
    [[maybe_unused]] const double dt) {
  assert(volFracs.size() == kNumberOfMaterials);
  assert(rhos.size() == kNumberOfMaterials);
  assert(ps.size() == kNumberOfMaterials);
  assert(soundSpeeds.size() == kNumberOfMaterials);
  assert(volFracRates.size() == kNumberOfMaterials);
  return 0.0;
}

#endif  // HYDRO_SOLVER_MATERIAL_CLOSURE_HPP_
