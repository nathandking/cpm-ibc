add_executable(DiffusionCurvesNefertiti DiffusionCurvesNefertiti.cpp)

target_link_libraries(DiffusionCurvesNefertiti
                        problems
                        src
                        )
                      
file(RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS DiffusionCurvesNefertiti RUNTIME DESTINATION ${REL})

set_target_properties(DiffusionCurvesNefertiti PROPERTIES FOLDER ${EXAMPLES_FOLDER})