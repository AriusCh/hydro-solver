#ifndef HYDRO_SOLVER_EOSES_MIGRUNEISEN_GENERAL_HPP_
#define HYDRO_SOLVER_EOSES_MIGRUNEISEN_GENERAL_HPP_

#include <cmath>

#include "../eos.hpp"

class EOSGruneisesenGeneral : public EOS {
 public:
  EOSGruneisesenGeneral(const double rho0_, const double C0_, const double G0_,
                        const double s_)
      : rho0(rho0_), C0(C0_), G0(G0_), s(s_) {}

 public:
  inline virtual double getp(double rho, double e) const override {
    const double eta = rho / rho0;
    const double output = rho0 * C0 * C0 * f(eta) + G0 * rho0 * e;
    return output;
  }
  inline virtual double gete(double rho, double p) const override {
    const double eta = rho / rho0;
    const double output = (p - rho0 * C0 * C0 * f(eta)) / (G0 * rho0);
    return output;
  }
  inline virtual double getc(double rho, double p) const override {
    const double eta = rho / rho0;
    const double output_squared =

        C0 * C0 * f_deta(eta) + G0 * p / (rho0 * eta * eta);
    return std::sqrt(output_squared);
  }

 private:
  double f(const double eta) const {
    const double output = (eta - 1) * (eta - 0.5 * G0 * (eta - 1)) /
                          std::pow((eta - s * (eta - 1)), 2);
    return output;
  }
  double f_deta(const double eta) const {
    const double output =
        (eta + (s - G0) * (eta - 1)) / std::pow(eta - s * (eta - 1), 3);
    return output;
  }

 private:
  const double rho0;
  const double C0;
  const double G0;
  const double s;
};

#endif  // HYDRO_SOLVER_EOSES_MIGRUNEISEN_GENERAL_HPP_
