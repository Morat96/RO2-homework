# Traveling Salesman Problem

## Requirements

CPLEX library and headers. <br>
Concorde library and headers.

## Build

mkdir build <br>
cd build <br>
cmake .. <br>
make

## Parameters
```
$ ./tsp -file FILE -time_limit TIME -model_type [0,7] -randomseed SEED -loop {0,1} -callback {0,1} -hardfixing {0,1} -localbranching {0,1} -vns {0,1} -ts {0,1} -sa {0,1} -ga {0,1}
```
- model_type
  - 0 -> standard model
  - 1 -> MTZ
  - 2 -> Flow 1
  - 3 -> Flow 2
  - 4 -> Flow 3
  - 5 -> T1
  - 6 -> T2
  - 7 -> T3

## Algorithms
- Compact Models
  - MTZ
  - Flow 1
  - Flow 2
  - Flow 3
  - T1
  - T2
  - T3

- Loop Methods
  - Simple Loop
  - Heuristic Loop

- Legacy Callbacks
  - Lazy Callback
  - UserCut Callback
  - Heuristic Callback

- Generic Callbacks
  - Lazy Callback
  - UserCut Callback
  - Heuristic Callback

- Matheuristics
  - Hard Fixing
  - Local Branching

- Heuristic solution builders
  - Nearest Neighbourhood
  - GRASP
  - Insertion
  - Insertion with Convex Hull

- Refining Algorithms
  - 2-OPT move
  - 3-OPT move
  
- Metaheuristics
  - Multi-start
  - VNS
  - Tabu Search
  - Simulated Annealing
  - Genetic Algorithm
