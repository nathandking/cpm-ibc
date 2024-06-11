add_executable(PoissonProblems PoissonProblems.cpp)

target_link_libraries(PoissonProblems
                        problems
                        src
                        )
                      
file(RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS PoissonProblems RUNTIME DESTINATION ${REL})

set_target_properties(PoissonProblems PROPERTIES FOLDER ${CONVERGENCE_STUDIES_FOLDER})