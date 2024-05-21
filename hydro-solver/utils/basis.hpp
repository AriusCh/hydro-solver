#ifndef HYDRO_SOLVER_UTILS_BASIS_HPP_
#define HYDRO_SOLVER_UTILS_BASIS_HPP_

#include <cstddef>

double lobattoBasis1D(const double x, const std::size_t order,
                      const std::size_t k);

double lobattoBasis1Ddx(const double x, const std::size_t order,
                        const std::size_t k);

double legendreBasis1D(const double x, const std::size_t order,
                       const std::size_t k);

double getBernsteinBasisCoeff(const std::size_t order, const std::size_t k);

double bernsteinBasis1D(const double x, const std::size_t order,
                        const std::size_t k);

double bernsteinBasis1Ddx(const double x, const std::size_t order,
                          const std::size_t k);

#endif  // HYDRO_SOLVER_UTILS_BASIS_HPP_
