add_executable(VectorFieldDesignSurfaceOfRevolution VectorFieldDesignSurfaceOfRevolution.cpp)

target_link_libraries(VectorFieldDesignSurfaceOfRevolution
                        problems
                        src
                        )
                      
file(RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS VectorFieldDesignSurfaceOfRevolution RUNTIME DESTINATION ${REL})

set_target_properties(VectorFieldDesignSurfaceOfRevolution PROPERTIES FOLDER ${EXAMPLES_FOLDER})