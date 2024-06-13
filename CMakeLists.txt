cmake_minimum_required(VERSION 3.25..3.29)

project(hydro-solver VERSION 0.1
  LANGUAGES CXX
  )

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_subdirectory(hydro-solver)
