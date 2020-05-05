# Ricerca Operativa 2
## Cartelle istanze
### Compact
Istanze per modelli compatti (devono avere pochi nodi, altrimenti ci mettono una vita, direi <= 25), si può eventualmente creare qualche istanza prendendo solo alcuni nodi da una istanza più grande. #instanze = almeno una decina.

### Dataset 
Istanze per algoritmi veloci ottimi (si può aggiungere anche qualche istanza più grande). #instanze = max 20.

### Heuristic
Istanze per algoritmi euristici (numero nodi >= 300). #instanze = max 20. <br>
Per ogni istanza -> 4/5 diversi randomseed. <br>
Ci dividiamo le run per tipo di algoritmo (es. compacts, heuristics ecc.).

## Parameters
```
$ ./tsp -file FILE -time_limit TIME -model_type [0,7] -randomseed SEED -loop {0,1} -callback {0,1} -hardfixing {0,1} -localbranching {0,1}
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
  - alone or with userCutCallback (Concorde)

- Generic callback
  - alone or with relaxation (Concorde)

- Hard Fixing (setting parameters)

- Local Branching (setting parameters)

- Heuristics
  - Nearest Neighbourhood (greedy)
  - Grasp (randomization)
  - Insertion
  - Insertion with Convex Hull




