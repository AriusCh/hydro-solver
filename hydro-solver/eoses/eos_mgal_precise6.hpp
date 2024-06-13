#ifndef HYDRO_SOLVER_EOSES_EOS_MGAL_PRECISE6_HPP_
#define HYDRO_SOLVER_EOSES_EOS_MGAL_PRECISE6_HPP_

#include <cmath>

#include "../eos.hpp"

class EOSMGAlPrecise6 : public EOS {
 public:
  EOSMGAlPrecise6()
      : rho0(2750.0), p0(560.964e9), a(1.12657), b(0.975511), ph(15e9) {}

 public:
  inline virtual double getp(double rho, double e) const override {
    const double x = rho / rho0;
    return pCold(rho) + G(x) * rho * (e - eCold(rho));
  }
  inline virtual double gete(double rho, double p) const override {
    const double x = rho / rho0;
    const double _ec = eCold(rho);
    const double _pc = pCold(rho);
    return _ec + 1. / rho / G(x) * (p - _pc);
  }
  inline virtual double getc(double rho, double p) const override {
    const double x = rho / rho0;
    const double _G = G(x);
    const double G_x = GPrime(x);
    const double pc = pCold(rho) / 1.e9;
    const double pcx = pColdPrime(rho) / 1.e9;
    const double c2 =
        1.e9 / rho0 *
        ((p / 1.e9 - pc) / (x * G(x)) * (_G + _G * _G + x * G_x) + pcx);

    return std::sqrt(c2);
  }

 private:
  double gets(double rho, double p) const;
  double getdpdrho(double rho, double e) const;
  double getdpde(double rho, double e) const;
  double G(double x) const;

  double GPrime(double x) const;
  double pCold(double rho) const;
  double pColdPrime(double rho) const;
  double eCold(double rho) const;
  double eColdPrime(double rho) const;

 private:
  const double rho0;
  const double p0;
  const double a;
  const double b;
  const double ph;
};

#endif  // HYDRO_SOLVER_EOSES_EOS_MGAL_PRECISE6_HPP_
