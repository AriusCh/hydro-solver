#include "pressure_relax_closure.hpp"

#include <cmath>

std::vector<double> PressureRelaxClosure::calcVolFracRates(
    const std::vector<double> &volFracs, const std::vector<double> &rhos,
    const std::vector<double> &ps, const std::vector<double> &soundSpeeds,
    [[maybe_unused]] const std::size_t cell, const double h, const double dt) {
  assert(volFracs.size() == kNumberOfMaterials);
  assert(rhos.size() == kNumberOfMaterials);
  assert(ps.size() == kNumberOfMaterials);
  assert(soundSpeeds.size() == kNumberOfMaterials);

  constexpr double Ctau = 0.25;
  constexpr double CL = 0.05;

  std::vector<double> bulkModuluses(kNumberOfMaterials, 0.0);
  for (std::size_t i = 0; i < kNumberOfMaterials; i++) {
    const double rho = rhos[i];
    const double soundSpeed = soundSpeeds[i];
    bulkModuluses[i] = rho * soundSpeed * soundSpeed;
  }

  double pStarNumerator = 0.0;
  double pStarDenominator = 0.0;
  for (std::size_t i = 0; i < kNumberOfMaterials; i++) {
    const double volFrac = volFracs[i];
    if (volFrac == 0) {
      continue;
    }
    const double p = ps[i];
    const double bulkModulus = bulkModuluses[i];
    const double div = volFrac / bulkModulus;
    pStarNumerator += p * div;
    pStarDenominator += div;
  }
  const double pStar = pStarNumerator / pStarDenominator;

  double minSoundSpeed = std::numeric_limits<double>::max();
  for (std::size_t i = 0; i < kNumberOfMaterials; i++) {
    const double volFrac = volFracs[i];
    if (volFrac == 0.0) {
      continue;
    }
    const double soundSpeed = soundSpeeds[i];
    if (soundSpeed < minSoundSpeed) {
      minSoundSpeed = soundSpeed;
    }
  }

  const double timeScale = Ctau * h / minSoundSpeed;

  std::vector<double> output(kNumberOfMaterials, 0.0);

  for (std::size_t i = 0; i < kNumberOfMaterials; i++) {
    const double volFrac = volFracs[i];
    if (volFrac == 0.0) {
      continue;
    }
    const double p = ps[i];
    const double bulkModulus = bulkModuluses[i];

    output[i] = (p - pStar) * volFrac / (bulkModulus * timeScale);
  }

  bool bShouldLimitRates = false;
  for (std::size_t i = 0; i < kNumberOfMaterials; i++) {
    const double volFrac = volFracs[i];
    if (volFrac == 0.0) {
      continue;
    }
    const double rate = output[i];
    const double absoluteChange = std::abs(dt * rate);
    const double maxChange = CL * volFrac;

    if (absoluteChange > maxChange) {
      bShouldLimitRates = true;
      break;
    }
  }
  if (bShouldLimitRates == false) {
    return output;
  }

  std::vector<double> maxRates(kNumberOfMaterials, 0.0);
  for (std::size_t i = 0; i < kNumberOfMaterials; i++) {
    const double volFrac = volFracs[i];
    if (volFrac == 0.0) {
      continue;
    }
    const double rate = output[i];
    const double maxRate = CL * volFrac / dt;
    const double maxRateSigned = std::copysign(maxRate, rate);
    maxRates[i] = maxRateSigned;
  }

  double ratePlus = 0.0;
  double rateMinus = 0.0;
  for (std::size_t i = 0; i < kNumberOfMaterials; i++) {
    const double volFrac = volFracs[i];
    if (volFrac == 0.0) {
      continue;
    }
    const double rate = maxRates[i];
    if (rate >= 0.0) {
      ratePlus += rate;
    } else {
      rateMinus += rate;
    }
  }

  for (std::size_t i = 0; i < kNumberOfMaterials; i++) {
    const double volFrac = volFracs[i];
    if (volFrac == 0.0) {
      continue;
    }
    const double rate = maxRates[i];
    if (rate > 0 && (ratePlus + rateMinus) > 0) {
      output[i] = -rateMinus / ratePlus * rate;
    } else if (rate < 0 && (ratePlus + rateMinus) < 0) {
      output[i] = -ratePlus / rateMinus * rate;
    } else {
      output[i] = rate;
    }
  }

  return output;
}

double PressureRelaxClosure::calcExchangePressure(
    const std::vector<double> &volFracs, const std::vector<double> &rhos,
    const std::vector<double> &ps,
    [[maybe_unused]] const std::vector<double> &soundSpeeds,
    const std::vector<double> &volFracRates, const double velocityScalarGrad,
    [[maybe_unused]] const double dt) {
  assert(volFracs.size() == kNumberOfMaterials);
  assert(rhos.size() == kNumberOfMaterials);
  assert(ps.size() == kNumberOfMaterials);
  assert(soundSpeeds.size() == kNumberOfMaterials);
  assert(volFracRates.size() == kNumberOfMaterials);

  double negativeRateMaxP = std::numeric_limits<double>::lowest();
  double positiveRateMinP = std::numeric_limits<double>::max();
  for (std::size_t material = 0; material < kNumberOfMaterials; material++) {
    const double volFrac = volFracs[material];
    if (volFrac == 0) {
      continue;
    }
    const double rate = volFracRates[material];
    const double p = ps[material];

    if (rate > 0 && p < positiveRateMinP) {
      positiveRateMinP = p;
    } else if (rate < 0 && p > negativeRateMaxP) {
      negativeRateMaxP = p;
    }
  }
  if (negativeRateMaxP == std::numeric_limits<double>::lowest() ||
      positiveRateMinP == std::numeric_limits<double>::max()) {
    return 0.0;
  }

  double negativeRateSum = 0.0;
  double positiveRateSum = 0.0;
  for (std::size_t material = 0; material < kNumberOfMaterials; material++) {
    const double volFrac = volFracs[material];
    if (volFrac == 0.0) {
      continue;
    }
    const double rate = volFracRates[material];
    const double rho = rhos[material];
    const double p = ps[material];

    negativeRateSum += rate * (p - negativeRateMaxP) / (volFrac * rho);
    positiveRateSum += rate * (p - positiveRateMinP) / (volFrac * rho);
  }

  const double maxSumP =
      positiveRateSum > negativeRateSum ? positiveRateMinP : negativeRateMaxP;
  const double minSumP =
      positiveRateSum > negativeRateSum ? negativeRateMaxP : positiveRateMinP;

  if (velocityScalarGrad <= 0.0) {
    return maxSumP;
  }
  return minSumP;
}
