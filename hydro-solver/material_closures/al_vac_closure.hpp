#ifndef HYDRO_SOLVER_MATERIAL_CLOSURES_AL_VAC_CLOSURE_HPP_
#define HYDRO_SOLVER_MATERIAL_CLOSURES_AL_VAC_CLOSURE_HPP_

#include "../material_closure.hpp"

class AlVacClosure : public MaterialClosure {
 public:
  AlVacClosure() : MaterialClosure(1) {}

  AlVacClosure(const AlVacClosure &) = delete;
  AlVacClosure(AlVacClosure &&) = delete;

  AlVacClosure &operator=(const AlVacClosure &) = delete;
  AlVacClosure &operator=(AlVacClosure &&) = delete;

  virtual ~AlVacClosure() = default;

 public:
  virtual std::vector<double> calcVolFracRates(
      const std::vector<double> &volFracs, const std::vector<double> &rhos,
      const std::vector<double> &ps, const std::vector<double> &soundSpeeds,
      const double h);
};

inline std::vector<double> AlVacClosure::calcVolFracRates(
    const std::vector<double> &volFracs, const std::vector<double> &rhos,
    const std::vector<double> &ps, const std::vector<double> &soundSpeeds,
    const double h) {
  assert(volFracs.size() == kNumberOfMaterials);
  assert(rhos.size() == kNumberOfMaterials);
  assert(ps.size() == kNumberOfMaterials);
  assert(soundSpeeds.size() == kNumberOfMaterials);

  constexpr double pSeparation = -1e9;
  constexpr double pStar = 0.0;
  constexpr double Ctau = 0.25;

  std::vector<double> output(kNumberOfMaterials, 0.0);

  const double pAl = ps[0];
  const double volFracAl = volFracs[0];

  if (volFracAl == 1.0 && pAl > pSeparation) {
    return output;
  }

  const double rhoAl = rhos[0];
  const double soundSpeedAl = soundSpeeds[0];

  const double timeScale = Ctau * h / soundSpeedAl;

  const double bulkModulus = rhoAl * soundSpeedAl * soundSpeedAl;

  const double rateAl =
      1.0 / timeScale * (pAl - pStar) * volFracAl / bulkModulus;

  output[0] = rateAl;

  return output;
}

#endif  // HYDRO_SOLVER_MATERIAL_CLOSURES_AL_VAC_CLOSURE_HPP_
