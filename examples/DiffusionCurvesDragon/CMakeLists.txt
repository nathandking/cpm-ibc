add_executable(DiffusionCurvesDragon DiffusionCurvesDragon.cpp)

target_link_libraries(DiffusionCurvesDragon
                        problems
                        src
                        )
                      
file(RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS DiffusionCurvesDragon RUNTIME DESTINATION ${REL})

set_target_properties(DiffusionCurvesDragon PROPERTIES FOLDER ${EXAMPLES_FOLDER})