#ifndef HYDRO_SOLVER_EOS_HPP_
#define HYDRO_SOLVER_EOS_HPP_

class EOS {
 public:
  EOS() = default;
  EOS(const EOS& rhs) = default;
  EOS(EOS&& rhs) = default;

  EOS& operator=(const EOS& rhs) = default;
  EOS& operator=(EOS&& rhs) = default;

  virtual ~EOS() = default;

 public:
  virtual double getp(const double rho, const double e) const = 0;
  virtual double gete(const double rho, const double p) const = 0;
  virtual double getc(const double rho, const double p) const = 0;
};

#endif  // HYDRO_SOLVER_EOS_HPP_
