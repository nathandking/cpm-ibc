add_executable(GeodesicDistanceDziukPoint GeodesicDistanceDziukPoint.cpp)

target_link_libraries(GeodesicDistanceDziukPoint
                        problems                
                        src
                        )
                      
file(RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS GeodesicDistanceDziukPoint RUNTIME DESTINATION ${REL})

set_target_properties(GeodesicDistanceDziukPoint PROPERTIES FOLDER ${EXAMPLES_FOLDER})