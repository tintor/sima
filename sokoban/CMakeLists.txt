project(sokoban)

add_executable(soko solver.cc corrals.cc level.cc main.cc)
target_link_libraries(soko core ${LIBRARIES})
