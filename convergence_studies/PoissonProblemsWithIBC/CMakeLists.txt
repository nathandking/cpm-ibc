add_executable(PoissonProblemsWithIBC PoissonProblemsWithIBC.cpp)

target_link_libraries(PoissonProblemsWithIBC
                        problems
                        src
                        )
                      
file(RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS PoissonProblemsWithIBC RUNTIME DESTINATION ${REL})

set_target_properties(PoissonProblemsWithIBC PROPERTIES FOLDER ${CONVERGENCE_STUDIES_FOLDER})