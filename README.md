# Traveling Salesman Problem

## Requirements

Link to CPLEX library and headers. <br>
Link to Concorde library and headers.

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
- Compacts
  - MTZ
  - Flow 1
  - Flow 2
  - Flow 3
  - T1
  - T2
  - T3

- Loop
  - I version
  - II version (two steps)

- Lazy callback
  - alone or with userCutCallback (Concorde) + heuristic callback + start MIP

- Generic callback
  - alone or with relaxation (Concorde) + heuristic callback + start MIP

- Hard Fixing (setting parameters)

- Local Branching (setting parameters)

- Heuristics
  - Nearest Neighbourhood (greedy)
  - Grasp (randomization)
  - Insertion
  - Insertion with Convex Hull

- Refining Algorithms
  - 2-OPT move
  - 3-OPT move
  
- Metaheuristics
  - VNS
  - Tabù Search
  - Simulated Annealing
  - Genetic Algorithm
  - Multi-start
  




