#ifndef HYDRO_SOLVER_MESH_GENERATOR_HPP_
#define HYDRO_SOLVER_MESH_GENERATOR_HPP_

#include <cstddef>

template <class Derived>
struct MeshGenerator {
  inline constexpr Derived& derived() { return *static_cast<Derived*>(this); }
  inline constexpr const Derived& derived() const {
    return *static_cast<Derived*>(this);
  }

  inline constexpr std::size_t GenerateMesh() {
    return derived().GenerateMesh();
  }
};

#endif  // HYDRO_SOLVER_MESH_GENERATOR_HPP_
