#ifndef HYDRO_SOLVER_MATERIAL_CLOSURES_AL_VAC_CELL_CLOSURE_HPP_
#define HYDRO_SOLVER_MATERIAL_CLOSURES_AL_VAC_CELL_CLOSURE_HPP_

#include <cmath>

#include "../material_closure.hpp"

class AlVacCellClosure : public MaterialClosure {
 public:
  AlVacCellClosure() : MaterialClosure(1) {}

  AlVacCellClosure(const AlVacCellClosure &) = delete;
  AlVacCellClosure(AlVacCellClosure &&) = delete;

  AlVacCellClosure &operator=(const AlVacCellClosure &) = delete;
  AlVacCellClosure &operator=(AlVacCellClosure &&) = delete;

  virtual ~AlVacCellClosure() = default;

 public:
  virtual std::vector<double> calcVolFracRates(
      const std::vector<double> &volFracs, const std::vector<double> &rhos,
      const std::vector<double> &ps, const std::vector<double> &soundSpeeds,
      const std::size_t cell, const double h, const double dt);

 private:
  enum class SeparationState { eWhole, eSeparated };
  std::vector<SeparationState> state = std::vector<SeparationState>(
      180000, AlVacCellClosure::SeparationState::eWhole);
};

inline std::vector<double> AlVacCellClosure::calcVolFracRates(
    const std::vector<double> &volFracs, const std::vector<double> &rhos,
    const std::vector<double> &ps, const std::vector<double> &soundSpeeds,
    const std::size_t cell, const double h, const double dt) {
  assert(volFracs.size() == kNumberOfMaterials);
  assert(rhos.size() == kNumberOfMaterials);
  assert(ps.size() == kNumberOfMaterials);
  assert(soundSpeeds.size() == kNumberOfMaterials);

  assert(cell < state.size());

  constexpr double pSeparation = -1e9;
  constexpr double pStar = 0.0;
  constexpr double Ctau = 1.0;
  constexpr double CL = 0.05;
  constexpr double DOUBLE_DELTA = 1e-14;

  std::vector<double> output(kNumberOfMaterials, 0.0);

  const double pAl = ps[0];
  const double volFracAl = volFracs[0];

  if (state[cell] == SeparationState::eWhole &&
      std::abs(volFracAl - 1.0) <= DOUBLE_DELTA && pAl > pSeparation) {
    return output;
  }
  state[cell] = SeparationState::eSeparated;

  const double rhoAl = rhos[0];
  const double soundSpeedAl = soundSpeeds[0];

  const double timeScale = Ctau * h / soundSpeedAl;

  const double bulkModulus = rhoAl * soundSpeedAl * soundSpeedAl;

  const double maxRateAl = CL * volFracAl / dt;

  double rateAl = 1.0 / timeScale * (pAl - pStar) * volFracAl / bulkModulus;

  if (volFracAl + rateAl * dt > 1.0 - DOUBLE_DELTA) {
    rateAl = (1.0 - DOUBLE_DELTA - volFracAl) / dt;
  }

  if (std::abs(rateAl) > maxRateAl) {
    rateAl = std::copysign(maxRateAl, rateAl);
  }

  output[0] = rateAl;

  return output;
}

#endif  // HYDRO_SOLVER_MATERIAL_CLOSURES_AL_VAC_CELL_CLOSURE_HPP_
