//
//  vns.h
//  cplex
//
//  Created by Matteo Moratello on 04/06/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//

#ifndef vns_h
#define vns_h

#include <stdio.h>
#include "tsp.h"

/**
 Variable Neighborhood Search Algorithm for TSP.
 
 @param inst instance of the struct "instance" for TSP problem.
 @param iter number of algorithm iterations.
 @param k value for the random K-move.
 */
void vns(instance* inst, int iter, int k);

//////////////////////////// utils ////////////////////////////

/**
 @brief Compute the distance between two nodes.
 
 @param i first node.
 @param j second node.
 @param inst instance of the struct "instance" for TSP problem.
 @return distance between two nodes.
 */
double dist(int i, int j, instance *inst);

/**
 Refining algorithm: 2-OPT move.
 
 @param inst instance of the struct "instance" for TSP problem.
 @param succ TPS solution with successors format.
 @param objval TSP objective function value.
 */
void twOpt(instance* inst, int *succ, double* objval);

/**
 Remove element from an array by sliding it.
 
 @param index array of int.
 @param from start index.
 @param to size of array.
 @return if removal was successful.
 */
int remove_index(int* index, int from, int to);

#endif /* vns_h */
