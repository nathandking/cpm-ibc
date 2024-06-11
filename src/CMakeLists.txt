add_library(src
            Tube.cpp
            IBCDOFSubset.cpp
            IBCMatrixManipulations.cpp
            IBCSubset.cpp
            IBCSubsets.cpp
            FDMatrices.cpp
            Geometry.cpp
            Grid.cpp
            Helpers.cpp
            Interpolation.cpp
            Node.cpp
            SimplePolygonMesh.cpp
            Surface.cpp
            SurfaceSpecifier.cpp
            Solver.cpp
            VectorMath.cpp
            Visualization.cpp)


target_link_libraries(src 
                      Eigen3::Eigen 
                      fcpw
                      lbfgspp)

if(USE_POLYSCOPE)
target_link_libraries(src polyscope)
target_compile_definitions(src PUBLIC POLYSCOPE) 
endif()

if(BUILD_ENABLE_SPARSE_GRID_SUPPORT)
    target_link_libraries(src 
                      nova::dynamics
                      nova::geometry
                      nova::spgrid
                      nova::tools
                    )
    target_compile_definitions(src PUBLIC SPGRID) 
endif()

if(NOT APPLE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
target_link_libraries(src OpenMP::OpenMP_CXX)
endif()

target_include_directories(src PUBLIC 
                            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                            $<INSTALL_INTERFACE:include>)

if(BUILD_CUSTOM_SOLVER)
    target_compile_definitions(src PUBLIC CUSTOM_SOLVER) 
else()
    if(BUILD_EIGEN_SPARSELU)
    target_compile_definitions(src PUBLIC ESL) 
    else()
    target_compile_definitions(src PUBLIC EBCG) 
    endif()
endif()