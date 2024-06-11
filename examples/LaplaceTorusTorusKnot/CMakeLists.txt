add_executable(LaplaceTorusTorusKnot LaplaceTorusTorusKnot.cpp)

target_link_libraries(LaplaceTorusTorusKnot
                        problems
                        src
                        )
                      
file(RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS LaplaceTorusTorusKnot RUNTIME DESTINATION ${REL})

set_target_properties(LaplaceTorusTorusKnot PROPERTIES FOLDER ${EXAMPLES_FOLDER})