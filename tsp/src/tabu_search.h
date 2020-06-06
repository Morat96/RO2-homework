//
//  tabu_search.h
//  cplex
//
//  Created by Matteo Moratello on 04/06/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//

#ifndef tabu_search_h
#define tabu_search_h

#include <stdio.h>
#include "tsp.h"

/**
 Structure that represent a couple of nodes, that is, an edge.
 */
typedef struct {
    int vert1, vert2;
} tabu_vert;

/**
 Tabù search algorithm for TSP.
 
 @param inst instance of the struct "instance" for TSP problem.
 @param iter number of iterations.
 @param list_size tabù list size.
 */
void tabu_search(instance* inst, int iter, int list_size);

/**
 Find the less worse edge swap in the TSP tour.
 
 @param inst instance of the struct "instance" for TSP problem.
 @param succ TPS solution with successors format.
 @param objval TSP objective function value.
 @param tabu_list tabù list.
 @param tabu_size current tabù list size.
 @param list_size original tabù list size.
 @param check if empty the list or not.
 */
void findBestEdgesSwap(instance* inst, int* succ, double* objval, tabu_vert* tabu_list, int* tabu_size, int list_size, int check);

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

#endif /* tabu_search_h */
