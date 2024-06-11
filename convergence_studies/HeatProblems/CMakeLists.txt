add_executable(HeatProblems HeatProblems.cpp)

target_link_libraries(HeatProblems
                        problems
                        src
                        )
                      
file(RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS HeatProblems RUNTIME DESTINATION ${REL})

set_target_properties(HeatProblems PROPERTIES FOLDER ${CONVERGENCE_STUDIES_FOLDER})