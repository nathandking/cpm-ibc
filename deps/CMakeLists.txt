set(DEPS_FOLDER deps)

add_subdirectory(fcpw)
if(USE_POLYSCOPE)
    add_subdirectory(polyscope)
endif()
add_subdirectory(nova)
add_subdirectory(LBFGSpp)
