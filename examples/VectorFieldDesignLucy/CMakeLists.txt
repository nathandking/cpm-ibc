add_executable(VectorFieldDesignLucy VectorFieldDesignLucy.cpp)

target_link_libraries(VectorFieldDesignLucy
                        problems
                        src
                        )
                      
file(RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS VectorFieldDesignLucy RUNTIME DESTINATION ${REL})

set_target_properties(VectorFieldDesignLucy PROPERTIES FOLDER ${EXAMPLES_FOLDER})