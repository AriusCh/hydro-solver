#include "basis.hpp"

#include <cmath>

#include "legendre_nodes.hpp"
#include "lobatto_nodes.hpp"

double lobattoBasis1D(const double x, const std::size_t order,
                      const std::size_t k) {
  assert(x >= 0.0 && x <= 1.0);
  assert(k <= order);

  const std::vector<double> nodes = getLobattoNodes(order);
  assert(nodes.size() == order + 1);

  double numerator{1.0};
  double denominator{1.0};

  for (std::size_t i = 0; i < nodes.size(); i++) {
    if (i == k) {
      continue;
    }

    numerator *= x - nodes[i];
    denominator *= nodes[k] - nodes[i];
  }

  const double output = numerator / denominator;
  return output;
}

double lobattoBasis1Ddx(const double x, const std::size_t order,
                        const std::size_t k) {
  assert(x >= 0.0 && x <= 1.0);
  assert(k <= order);

  std::vector<double> nodes = getLobattoNodes(order);

  double output{0.0};
  for (std::size_t i = 0; i < nodes.size(); i++) {
    if (i == k) {
      continue;
    }

    double numerator{1.0};
    double denominator{1.0};
    for (std::size_t j = 0; j < nodes.size(); j++) {
      if (j == i || j == k) {
        continue;
      }

      numerator *= x - nodes[j];
      denominator *= nodes[k] - nodes[j];
    }
    output += numerator / (denominator * (nodes[k] - nodes[i]));
  }

  return output;
}

double legendreBasis1D(const double x, const std::size_t order,
                       const std::size_t k) {
  assert(x >= 0.0 && x <= 1.0);
  assert(k <= order);

  const std::vector<double> nodes = getLegendreNodes(order);

  double numerator{1.0};
  double denominator{1.0};

  for (std::size_t i = 0; i < nodes.size(); i++) {
    if (i == k) {
      continue;
    }

    numerator *= x - nodes[i];
    denominator *= nodes[k] - nodes[i];
  }

  const double output = numerator / denominator;
  return output;
}

double getBernsteinBasisCoeff(const std::size_t order, const std::size_t k) {
  assert(k <= order);

  const std::size_t kmin = std::min(k, order - k);
  double numerator{1.0};
  double denominator{1.0};
  for (std::size_t i = 1; i <= kmin; i++) {
    numerator *= order - kmin + i;
    denominator *= i;
  }

  return numerator / denominator;
}

double bernsteinBasis1D(const double x, const std::size_t order,
                        const std::size_t k) {
  assert(x >= 0.0 && x <= 1.0);
  assert(k <= order);

  const double coeff = getBernsteinBasisCoeff(order, k);
  return coeff * std::pow(x, k) * std::pow(1.0 - x, order - k);
}

double bernsteinBasis1Ddx(const double x, const std::size_t order,
                          const std::size_t k) {
  assert(x >= 0.0 && x <= 1.0);
  assert(k <= order);

  const double coeff = getBernsteinBasisCoeff(order, k);

  double leftPart =
      k > 0 ? k * std::pow(x, k - 1) * std::pow(1.0 - x, order - k) : 0.0;
  double rightPart = k < order ? (order - k) * std::pow(x, k) *
                                     std::pow(1.0 - x, order - k - 1)
                               : 0.0;
  return coeff * (leftPart - rightPart);
}
