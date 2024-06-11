cmake_minimum_required(VERSION 3.11)

project(cpm-base VERSION 0.1.0 LANGUAGES CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
 
find_package(Eigen3 3.3 REQUIRED NO_MODULE) 
if(NOT APPLE)
find_package(OpenMP REQUIRED)
endif()

SET(FCPW_USE_ENOKI OFF CACHE BOOL "Build enoki")
SET(USE_POLYSCOPE ON CACHE BOOL "Build Polyscope")

add_subdirectory(convergence_studies)
add_subdirectory(deps)
add_subdirectory(examples)
add_subdirectory(problems)
add_subdirectory(src)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)   


set(BUILD_CUSTOM_SOLVER TRUE CACHE BOOL "use matrix-free solver")
if(NOT BUILD_CUSTOM_SOLVER)
set(BUILD_EIGEN_SPARSELU FALSE CACHE BOOL "use Eigen::SparseLU; otherwise Eigen::BiCGStab will be used")
endif()