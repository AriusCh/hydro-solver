#ifndef HYDRO_SOLVER_EOSES_EOS_LOGARITHMIC_HPP_
#define HYDRO_SOLVER_EOSES_EOS_LOGARITHMIC_HPP_

#include <cmath>

#include "../eos.hpp"

class EOSLogarithmic : public EOS {
 public:
  EOSLogarithmic() {}
  EOSLogarithmic(const EOSLogarithmic& rhs) = default;
  EOSLogarithmic(EOSLogarithmic&& rhs) = default;

  EOSLogarithmic& operator=(const EOSLogarithmic& rhs) = default;
  EOSLogarithmic& operator=(EOSLogarithmic&& rhs) = default;

  virtual ~EOSLogarithmic() = default;

 public:
  inline virtual double getp(double rho,
                             [[maybe_unused]] double e) const override {
    return B * std::log(rho / rho0);
  }
  inline virtual double gete([[maybe_unused]] double rho,
                             [[maybe_unused]] double p) const override {
    return 0.0;
  }
  inline virtual double getc(double rho,
                             [[maybe_unused]] double p) const override {
    return std::sqrt(B / rho);
  }

 private:
  static constexpr double rho0 = 2688.9;
  static constexpr double B = 73e9;
};

#endif  // HYDRO_SOLVER_EOSES_EOS_LOGARITHMIC_HPP_
