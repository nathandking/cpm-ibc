add_executable(GeodesicDistanceBeetle GeodesicDistanceBeetle.cpp)

target_link_libraries(GeodesicDistanceBeetle
                        problems                
                        src
                        )
                      
file(RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS GeodesicDistanceBeetle RUNTIME DESTINATION ${REL})

set_target_properties(GeodesicDistanceBeetle PROPERTIES FOLDER ${EXAMPLES_FOLDER})