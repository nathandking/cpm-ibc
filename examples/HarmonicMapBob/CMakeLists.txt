add_executable(HarmonicMapBob HarmonicMapBob.cpp)

target_link_libraries(HarmonicMapBob
                        problems
                        src
                        )
                      
file(RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS HarmonicMapBob RUNTIME DESTINATION ${REL})

set_target_properties(HarmonicMapBob PROPERTIES FOLDER ${EXAMPLES_FOLDER})