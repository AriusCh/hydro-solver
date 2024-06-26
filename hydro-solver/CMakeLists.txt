project(hydro-solver)

set(SOURCE_FILES main.cpp

  # EOSes
  eos.hpp
  eoses/eos_ideal_gas.hpp
  eoses/eos_mgal_precise6.hpp eoses/eos_mgal_precise6.cpp

  # MATERIAL CLOSURES
  material_closure.hpp
  material_closures/al_vac_closure.hpp

  # PROBLEM
  problem.hpp problem.cpp

  # METHODS
  method.hpp
  methods/lagrangian_fem_method.hpp methods/lagrangian_fem_method.cpp

  # SIMULATION
  simulation.hpp simulation.cpp

  # UTILS
  utils/lobatto_nodes.hpp utils/lobatto_nodes.cpp
  utils/legendre_nodes.hpp utils/legendre_nodes.cpp
  utils/legendre_weights.hpp utils/legendre_weights.cpp
  utils/basis.hpp utils/basis.cpp
  )

add_executable(hydro-solver ${SOURCE_FILES})

set_target_properties(hydro-solver PROPERTIES
  CXX_STANDARD 17
  CXX_STANDARD_REQUIRED YES
  CXX_EXTENSIONS NO
  )

find_package(Eigen3 3.4 REQUIRED NO_MODULE)
target_link_libraries(hydro-solver PUBLIC Eigen3::Eigen)

find_package(OpenMP)
if (OpenMP_CXX_FOUND)
  target_link_libraries(hydro-solver PUBLIC OpenMP::OpenMP_CXX)
endif()

target_compile_options(hydro-solver PRIVATE -Wall -Wextra -Wpedantic -Wshadow)
