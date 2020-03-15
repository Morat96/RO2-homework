cmake_minimum_required (VERSION 2.8)

PROJECT (TSP)

FIND_PACKAGE(CPLEX REQUIRED)

INCLUDE_DIRECTORIES(include ${CPLEX_INCLUDE_DIRS})

set(SRC
    src/tsp.c
    src/tsp.h
    src/main.c
)

set (CMAKE_CXX_STANDARD 11)

add_executable(tsp ${SRC})

target_link_libraries(tsp ${CPLEX_LIBRARY})

set(CMAKE_SUPPRESS_REGENERATION true)