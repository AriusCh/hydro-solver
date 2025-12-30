#ifndef HYDRO_SOLVER_EOSES_EOS_SIMPLE_PLASTIC_HPP_
#define HYDRO_SOLVER_EOSES_EOS_SIMPLE_PLASTIC_HPP_

#include <cmath>

#include "../eos.hpp"

class EOSSimplePlastic : public EOS {
 public:
  EOSSimplePlastic(const double bulk_modulus_, const double rho0_)
      : bulk_modulus(bulk_modulus_), rho0(rho0_) {}
  EOSSimplePlastic(const EOSSimplePlastic& rhs) = default;
  EOSSimplePlastic(EOSSimplePlastic&& rhs) = default;

  EOSSimplePlastic& operator=(const EOSSimplePlastic& rhs) = default;
  EOSSimplePlastic& operator=(EOSSimplePlastic&& rhs) = default;

  virtual ~EOSSimplePlastic() = default;

 public:
  inline virtual double getp(double rho,
                             [[maybe_unused]] double e) const override {
    return bulk_modulus * (rho - rho0) / rho0;
  }
  inline virtual double gete([[maybe_unused]] double rho,
                             [[maybe_unused]] double p) const override {
    return 0.0;
  }
  inline virtual double getc([[maybe_unused]] double rho,
                             [[maybe_unused]] double p) const override {
    return std::sqrt(bulk_modulus / rho0);
  }

 private:
  double bulk_modulus;
  double rho0;
};

#endif  // HYDRO_SOLVER_EOSES_EOS_SIMPLE_PLASTIC_HPP_
