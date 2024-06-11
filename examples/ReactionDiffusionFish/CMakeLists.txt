add_executable(ReactionDiffusionFish ReactionDiffusionFish.cpp)

target_link_libraries(ReactionDiffusionFish
                        problems
                        src
                        )
                      
file(RELATIVE_PATH REL ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS ReactionDiffusionFish RUNTIME DESTINATION ${REL})

set_target_properties(ReactionDiffusionFish PROPERTIES FOLDER ${EXAMPLES_FOLDER})