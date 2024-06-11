add_executable(HeatProblemsWithIBC HeatProblemsWithIBC.cpp)

target_link_libraries(HeatProblemsWithIBC
                        problems
                        src
                        )
                      
file(RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS HeatProblemsWithIBC RUNTIME DESTINATION ${REL})

set_target_properties(HeatProblemsWithIBC PROPERTIES FOLDER ${CONVERGENCE_STUDIES_FOLDER})