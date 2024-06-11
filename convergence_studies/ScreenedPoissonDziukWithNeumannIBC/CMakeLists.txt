add_executable(ScreenedPoissonDziukWithNeumannIBC ScreenedPoissonDziukWithNeumannIBC.cpp)

target_link_libraries(ScreenedPoissonDziukWithNeumannIBC
                        problems                
                        src
                        )
                      
file(RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS ScreenedPoissonDziukWithNeumannIBC RUNTIME DESTINATION ${REL})

set_target_properties(ScreenedPoissonDziukWithNeumannIBC PROPERTIES FOLDER ${CONVERGENCE_STUDIES_FOLDER})