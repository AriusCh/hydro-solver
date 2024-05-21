#include "legendre_weights.hpp"

#include <stdexcept>
#include <string>

std::vector<double> getLegendreWeights(const std::size_t order) {
  if (order > legendreOrderMax) {
    const std::string error_message{
        "REQUESTED LEGENDRE ORDER: " + std::to_string(order) +
        ", MAX LEGENDRE ORDER: " + std::to_string(legendreOrderMax)};

    throw std::runtime_error(error_message);
  }

  const std::size_t indMin = getLegendreStartIndex(order);
  std::vector<double> output{legendreWeights.cbegin() + indMin,
                             legendreWeights.cbegin() + indMin + order + 1};

  return output;
}
