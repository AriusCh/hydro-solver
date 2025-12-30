#ifndef HYDRO_SOLVER_EOSES_EOS_MGAL_SOLID_HPP_
#define HYDRO_SOLVER_EOSES_EOS_MGAL_SOLID_HPP_

#include <cmath>
#include <iostream>

#include "../eos.hpp"

class EOSMGAlSolid : public EOS {
 public:
  EOSMGAlSolid() {}
  // : rho0(2730.0), a(1.12657), b(0.975511), p0(560.964e9) {}

 public:
  inline virtual double getp(double rho, double e) const override {
    double x = rho / rho0;
    return p0 * (pCold(rho) + G(x) * x * (e / e0 - eCold(rho)));
  }

  inline virtual double gete(double rho, double p) const override {
    double x = rho / rho0;
    double _ec = eCold(rho);
    double _pc = pCold(rho);
    return ((p / p0 - _pc) / G(x) / x + _ec) * e0;
  }

  inline virtual double getc(double rho, double p) const override {
    constexpr double eps = 1e-10;
    const double x = rho / rho0;
    const double e = gete(rho, p);
    const double dpdr =
        (getp(rho * (1.0 + eps), e) - getp(rho, e)) / (rho * eps);
    const double dpde = G(x) * x / e0 * p0;
    const double c_sqr = dpdr + p * dpde / rho / rho;

    if (c_sqr < 10.0) {
      std::cout << "Fuck\n";
      std::cout << "rho: " << rho;
      std::cout << "   p: " << p;
      std::cout << "   e: " << e << '\n';
      std::cout << "dpdr: " << dpdr;
      std::cout << "   dpde: " << dpde;
      std::cout << "   c_sqr: " << c_sqr << '\n';

      exit(EXIT_FAILURE);
    }

    return std::sqrt(c_sqr);
  }

 private:
  inline double G([[maybe_unused]] const double x) const { return 0.5; }

  double pCold(const double rho) const {
    const double x = rho / rho0;
    const double pc = p_split * 0.5 * (1.0 - std::copysign(1.0, x - x_split)) +
                      x * (std::pow(x, a) - std::pow(x, b)) * 0.5 *
                          (1.0 + std::copysign(1.0, x - x_split));

    return pc;
  }
  double eCold(const double rho) const {
    const double x = rho / rho0;
    const double ec = (e_split + p_split * (1.0 / x_split - 1.0 / x)) * 0.5 *
                          (1.0 - std::copysign(1.0, x - x_split)) +
                      (std::pow(x, a) / a - std::pow(x, b) / b) * 0.5 *
                          (1.0 + std::copysign(1.0, x - x_split));
    return ec;
  }

 private:
  static constexpr double rho0 = 2730.0;
  static constexpr double a = 1.12657;
  static constexpr double b = 0.975511;
  static constexpr double p0 = 560.964e9;
  static constexpr double e0 = p0 / rho0;

  static constexpr double x_split = 0.613989;
  static constexpr double p_split = -0.027100423868395146;
  static constexpr double e_split = -0.12458854272720732;
};

#endif  // HYDRO_SOLVER_EOSES_EOS_MGAL_SOLID_HPP_
