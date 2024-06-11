add_executable(DiffusionCurvesCodimZero DiffusionCurvesCodimZero.cpp)

target_link_libraries(DiffusionCurvesCodimZero
                        problems
                        src
                        )
                      
file(RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS DiffusionCurvesCodimZero RUNTIME DESTINATION ${REL})

set_target_properties(DiffusionCurvesCodimZero PROPERTIES FOLDER ${EXAMPLES_FOLDER})