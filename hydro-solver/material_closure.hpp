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
      const double h);

 protected:
  const std::size_t kNumberOfMaterials;
};

inline std::vector<double> MaterialClosure::calcVolFracRates(
    [[maybe_unused]] const std::vector<double> &volFracs,
    [[maybe_unused]] const std::vector<double> &rhos,
    [[maybe_unused]] const std::vector<double> &ps,
    [[maybe_unused]] const std::vector<double> &soundSpeeds,
    [[maybe_unused]] const double h) {
  assert(volFracs.size() == kNumberOfMaterials);
  assert(rhos.size() == kNumberOfMaterials);
  assert(ps.size() == kNumberOfMaterials);
  assert(soundSpeeds.size() == kNumberOfMaterials);
  return std::vector<double>(kNumberOfMaterials, 0.0);
}

#endif  // HYDRO_SOLVER_MATERIAL_CLOSURE_HPP_
