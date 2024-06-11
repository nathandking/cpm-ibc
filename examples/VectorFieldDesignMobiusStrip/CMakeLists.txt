add_executable(VectorFieldDesignMobiusStrip VectorFieldDesignMobiusStrip.cpp)

target_link_libraries(VectorFieldDesignMobiusStrip
                        problems
                        src
                        )
                      
file(RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS VectorFieldDesignMobiusStrip RUNTIME DESTINATION ${REL})

set_target_properties(VectorFieldDesignMobiusStrip PROPERTIES FOLDER ${EXAMPLES_FOLDER})