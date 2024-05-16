project(hydro-solver)

set(SOURCE_FILES main.cpp)

add_executable(hydro-solver ${SOURCE_FILES})

set_target_properties(hydro-solver PROPERTIES
  CXX_STANDARD 23
  CXX_STANDARD_REQUIRED YES
  CXX_EXTENSIONS NO
  )

find_package(Eigen3 3.4 REQUIRED NO_MODULE)
target_link_libraries(hydro-solver PUBLIC Eigen3::Eigen)

# find_package(OpenMP)
# if (OpenMP_CXX_FOUND)
#   target_link_libraries(hydro-solver PUBLIC OpenMP::OpenMP_CXX)
# endif()

target_compile_options(hydro-solver PRIVATE -Wall -Wextra -Wpedantic -Wshadow -Weffc++)
