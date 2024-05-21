#ifndef HYDRO_SOLVER_EOSES_EOS_IDEAL_GAS_HPP_
#define HYDRO_SOLVER_EOSES_EOS_IDEAL_GAS_HPP_

#include <cmath>

#include "../eos.hpp"

class EOSIdealGas : public EOS {
 public:
  EOSIdealGas(const double gamma_) : gamma(gamma_) {}
  EOSIdealGas(const EOSIdealGas& rhs) = default;
  EOSIdealGas(EOSIdealGas&& rhs) = default;

  EOSIdealGas& operator=(const EOSIdealGas& rhs) = default;
  EOSIdealGas& operator=(EOSIdealGas&& rhs) = default;

  virtual ~EOSIdealGas() = default;

 public:
  inline virtual double getp(double rho, double e) const override {
    return (gamma - 1.0) * rho * e;
  }
  inline virtual double gete(double rho, double p) const override {
    return p / ((gamma - 1.0) * rho);
  }
  inline virtual double getc(double rho, double p) const override {
    return std::sqrt(gamma * p / rho);
  }

 private:
  double gamma;
};

class EOSIdealGasWithBackgroundPressure : public EOS {
 public:
  EOSIdealGasWithBackgroundPressure(const double gamma_, const double p0_)
      : gamma(gamma_), p0(p0_) {}
  EOSIdealGasWithBackgroundPressure(
      const EOSIdealGasWithBackgroundPressure& rhs) = default;
  EOSIdealGasWithBackgroundPressure(EOSIdealGasWithBackgroundPressure&& rhs) =
      default;

  EOSIdealGasWithBackgroundPressure& operator=(
      const EOSIdealGasWithBackgroundPressure& rhs) = default;
  EOSIdealGasWithBackgroundPressure& operator=(
      EOSIdealGasWithBackgroundPressure&& rhs) = default;

  virtual ~EOSIdealGasWithBackgroundPressure() = default;

 public:
  inline virtual double getp(const double rho, const double e) const override {
    return (gamma - 1.0) * rho * e - gamma * p0;
  }
  inline virtual double gete(const double rho, const double p) const override {
    return (p + gamma * p0) / ((gamma - 1.0) * rho);
  }
  inline virtual double getc(const double rho, const double p) const override {
    return std::sqrt(gamma * std::abs(p + gamma * p0) / rho);
  }

 private:
  double gamma;
  double p0;
};

#endif  // HYDRO_SOLVER_EOSES_EOS_IDEAL_GAS_HPP_
