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
    src/hardfixing.h
    src/hardfixing.c
    src/localbranching.h
    src/localbranching.c
    src/loop_method.h
    src/loop_method.c
    src/compact_models.h
    src/compact_models.c
)

set (CMAKE_CXX_STANDARD 11)

add_executable(tsp ${SRC})

target_link_libraries(tsp ${CPLEX_LIBRARY})
target_link_libraries(tsp ${CONCORDE_LIBRARY})

set(CMAKE_SUPPRESS_REGENERATION true)