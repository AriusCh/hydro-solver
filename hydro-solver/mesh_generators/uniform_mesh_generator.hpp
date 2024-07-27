#ifndef HYDRO_SOLVER_UNIFORM_MESH_GENERATOR_HPP_
#define HYDRO_SOLVER_UNIFORM_MESH_GENERATOR_HPP_

#include <array>
#include <cstddef>

#include "../mesh_generator.hpp"

template <std::size_t dimension>
class UniformMeshGenerator
    : public MeshGenerator<UniformMeshGenerator<dimension>> {
 public:
  inline constexpr std::size_t GenerateMesh() { return 0; }

 public:
  std::array<double, dimension> mins;
  std::array<double, dimension> maxes;
};

#endif  // HYDRO_SOLVER_UNIFORM_MESH_GENERATOR_HPP_
