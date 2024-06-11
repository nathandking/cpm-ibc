if(NOT APPLE)
FIND_PACKAGE(OpenMP)
  if(OPENMP_FOUND)
     message("OPENMP FOUND")
     set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
     set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif()
 include_directories(SYSTEM ${OpenMP_INCLUDE_PATH})
endif()

add_library(problems
            GeodesicDistance.cpp
            GeodesicDistanceNearest.cpp
            Heat.cpp
            HeatWithIBC.cpp
            Poisson.cpp
            PoissonWithIBC.cpp
            PoissonWithIBCNearest.cpp
            ReactionDiffusionWithIBC.cpp)

target_link_libraries(problems 
                      src
                      Eigen3::Eigen 
                      fcpw
                     )

if(USE_POLYSCOPE)
target_link_libraries(problems polyscope)
target_compile_definitions(problems PUBLIC POLYSCOPE) 
endif()

target_include_directories(problems PUBLIC 
                            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                            $<INSTALL_INTERFACE:include>)

if(BUILD_CUSTOM_SOLVER)
    target_compile_definitions(problems PUBLIC CUSTOM_SOLVER) 
else()
    if(BUILD_EIGEN_SPARSELU)
    target_compile_definitions(problems PUBLIC ESL) 
    else()
    target_compile_definitions(problems PUBLIC EBCG) 
    endif()
endif()