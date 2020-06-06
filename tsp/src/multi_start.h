//
//  multi_start.h
//  cplex
//
//  Created by Matteo Moratello on 04/06/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//

#ifndef multi_start_h
#define multi_start_h

#include <stdio.h>
#include "tsp.h"

/**
 Multi-start algorithm for TSP.
 
 @param inst instance of the struct "instance" for TSP problem.
 @param iter number of iterations.
 */
void multi_start(instance* inst, int iter);

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
 Build a random TSP feasible solution.
 
 @param inst instance of the struct "instance" for TSP problem.
 @param succ TSP tour represented with successors.
 @param objval TSP objective function value.
 */
void random_solution_ms(instance* inst, int* succ, double* objval);

/**
 Build a random TSP feasible solution using Grasp.
 
 @param inst inst instance of the struct "instance" for TSP problem.
 @param succ succ TSP tour represented with successors.
 @param objval objval TSP objective function value.
 */
void grasp_solution_ms(instance *inst, int* succ, double* objval);

/**
 Refining algorithm: 2-OPT move.
 
 @param inst instance of the struct "instance" for TSP problem.
 @param succ TPS solution with successors format.
 @param objval TSP objective function value.
 */
void twOpt(instance* inst, int *succ, double* objval);

/**
 Compute minimum of an array and return the corresponding index.
 
 @param array array of double.
 @param arr_size array size.
 @return array index of minimum value of the array.
 */
int min(double *array, int arr_size);

/**
 Compute the three smaller values of an array and store the indexes in the "ind" array.
 
 @param array array of double.
 @param arr_size array size.
 @param ind array of int containing the three smallest values of array.
 @return return 1 if the array size is smaller than three, otherwise return 0.
 */
int three_min(double *array, int arr_size, int *ind);

#endif /* multi_start_h */
