#include "lobatto_nodes.hpp"

#include <stdexcept>
#include <string>

std::vector<double> getLobattoNodes(const std::size_t order) {
  if (order > lobattoOrderMax) {
    const std::string error_message{
        "REQUESTED LOBATTO ORDER: " + std::to_string(order) +
        ", MAX LOBATTO ORDER: " + std::to_string(lobattoOrderMax)};

    throw std::runtime_error(error_message);
  }

  const std::size_t indMin = getLobattoStartIndex(order);
  std::vector<double> output{lobattoNodes.cbegin() + indMin,
                             lobattoNodes.cbegin() + indMin + order + 1};

  return output;
}
