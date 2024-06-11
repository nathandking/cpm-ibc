add_executable(DiffusionCurvesMixedCodim DiffusionCurvesMixedCodim.cpp)

target_link_libraries(DiffusionCurvesMixedCodim
                        problems
                        src
                        )
                      
file(RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS DiffusionCurvesMixedCodim RUNTIME DESTINATION ${REL})

set_target_properties(DiffusionCurvesMixedCodim PROPERTIES FOLDER ${EXAMPLES_FOLDER})