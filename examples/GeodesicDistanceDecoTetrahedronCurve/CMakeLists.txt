add_executable(GeodesicDistanceDecoTetrahedronCurve GeodesicDistanceDecoTetrahedronCurve.cpp)

target_link_libraries(GeodesicDistanceDecoTetrahedronCurve
                        problems                
                        src
                        )
                      
file(RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS GeodesicDistanceDecoTetrahedronCurve RUNTIME DESTINATION ${REL})

set_target_properties(GeodesicDistanceDecoTetrahedronCurve PROPERTIES FOLDER ${EXAMPLES_FOLDER})