add_executable(GeodesicDistanceTorusTrefoil GeodesicDistanceTorusTrefoil.cpp)

target_link_libraries(GeodesicDistanceTorusTrefoil
                        problems                
                        src
                        )
                      
file(RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS GeodesicDistanceTorusTrefoil RUNTIME DESTINATION ${REL})

set_target_properties(GeodesicDistanceTorusTrefoil PROPERTIES FOLDER ${EXAMPLES_FOLDER})