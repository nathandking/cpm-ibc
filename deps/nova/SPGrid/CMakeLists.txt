set(BUILD_ENABLE_SPARSE_GRID_SUPPORT OFF CACHE BOOL "Enable building the SPGrid framework.")

if(BUILD_ENABLE_SPARSE_GRID_SUPPORT)
# Target for SPGrid
add_library( _nova_spgrid
  src/Core/SPGrid_Geometry.cpp
  src/Core/SPGrid_Utilities.cpp
  )

target_include_directories( _nova_spgrid PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../Geometry/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE src )

# Add required compiler features
# See here for more information:
# https://cmake.org/cmake/help/v3.1/prop_gbl/CMAKE_CXX_KNOWN_FEATURES.html#prop_gbl:CMAKE_CXX_KNOWN_FEATURES

target_compile_features( _nova_spgrid
  PRIVATE
  cxx_nullptr
  )

set(SPARSE_GRID_USE_HASWELL FALSE CACHE BOOL "use haswell-enabled bit-interleaving hardware instruction.")

if(SPARSE_GRID_USE_HASWELL)
  if("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    target_compile_options( _nova_spgrid PUBLIC -mbmi2)
  endif()
  target_compile_definitions( _nova_spgrid PUBLIC HASWELL=1)
endif()

set(USE_SPARSE_GRID FALSE CACHE BOOL "use sparse grid")
if(USE_SPARSE_GRID)
target_compile_definitions(_nova_spgrid PUBLIC  
USE_SPARSE_GRID  
) 
endif()

target_link_libraries( _nova_spgrid
  PUBLIC
  nova::tools
  )


add_library( nova::spgrid ALIAS _nova_spgrid )
endif()
