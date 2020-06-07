cmake_minimum_required (VERSION 2.8)

PROJECT (TSP)

FIND_PACKAGE(CPLEX REQUIRED)

INCLUDE_DIRECTORIES(include ${CPLEX_INCLUDE_DIRS})

FIND_PATH(CONCORDE_INCLUDE_DIR
	cut.h
  	PATHS "/Users/matteomoratello/Downloads/concorde/INCLUDE/")

FIND_LIBRARY(CONCORDE_LIBRARY
  	concorde.a
  	PATHS "/Users/matteomoratello/Downloads/concorde/")

INCLUDE_DIRECTORIES(include ${CONCORDE_INCLUDE_DIR})

set(SRC
    src/tsp.c
    src/tsp.h
    src/utilities.h
    src/utilities.c
    src/main.c
    src/chrono.c
    src/callback.h
    src/callback.c
    src/convexhull.cpp
    src/heuristics.h
    src/heuristics.c
    src/hardfixing.h
    src/hardfixing.c
    src/localbranching.h
    src/localbranching.c
    src/loop_method.h
    src/loop_method.c
    src/compact_models.h
    src/compact_models.c
    src/vns.h
    src/vns.c
    src/tabu_search.h
    src/tabu_search.c
    src/simulated_annealing.h
    src/simulated_annealing.c
    src/genetic_algorithm.h
    src/genetic_algorithm.c
    src/multi_start.h
    src/multi_start.c
)

set (CMAKE_CXX_STANDARD 11)

add_executable(tsp ${SRC})

target_link_libraries(tsp ${CPLEX_LIBRARY})
target_link_libraries(tsp ${CONCORDE_LIBRARY})

set(CMAKE_SUPPRESS_REGENERATION true)