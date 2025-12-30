#ifndef HYDRO_SOLVER_EOSES_STIFFENED_GAS_HPP_
#define HYDRO_SOLVER_EOSES_STIFFENED_GAS_HPP_

#include <cmath>

#include "../eos.hpp"

class EOSStiffenedGas : public EOS {
 public:
  EOSStiffenedGas(const double rho0_, const double c0_, const double g_,
                  const double bulk_modulus_)
      : rho0(rho0_), c0(c0_), g(g_), bulk_modulus(bulk_modulus_) {}

 public:
  inline virtual double getp(double rho, double e) const override {
    const double p = c0 * c0 * (rho - rho0) + (g - 1) * rho * e;
    return p;
  }

  inline virtual double gete(double rho, double p) const override {
    const double e = (p - c0 * c0 * (rho - rho0)) / ((g - 1) * rho);
    return e;
  }

  inline virtual double getc(double rho,
                             [[maybe_unused]] double p) const override {
    const double c = std::sqrt(bulk_modulus / rho);
    return c;
  }

 private:
  const double rho0;
  const double c0;
  const double g;
  const double bulk_modulus;
};

#endif  // HYDRO_SOLVER_EOSES_STIFFENED_GAS_HPP_
