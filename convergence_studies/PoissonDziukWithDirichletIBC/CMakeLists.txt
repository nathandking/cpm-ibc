add_executable(PoissonDziukWithDirichletIBC PoissonDziukWithDirichletIBC.cpp)

target_link_libraries(PoissonDziukWithDirichletIBC
                        problems                
                        src
                        )
                      
file(RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS PoissonDziukWithDirichletIBC RUNTIME DESTINATION ${REL})

set_target_properties(PoissonDziukWithDirichletIBC PROPERTIES FOLDER ${CONVERGENCE_STUDIES_FOLDER})