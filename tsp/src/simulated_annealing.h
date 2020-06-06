//
//  simulated_annealing.h
//  cplex
//
//  Created by Matteo Moratello on 04/06/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//

#ifndef simulated_annealing_h
#define simulated_annealing_h

#include <stdio.h>
#include "tsp.h"

/**
 Simulated annealing algorithm for TSP.
 
 @param inst instance of the struct "instance" for TSP problem.
 @param iter number of iterations.
 @param size temperature list size.
 */
void simulated_annealing(instance* inst, int iter, int size);

/**
 Swap two random edges of TSP tour.
 
 @param inst instance of the struct "instance" for TSP problem.
 @param succ TPS solution with successors format.
 @param objval TSP objective function value.
 */
void neigh_random_swap(instance* inst, int* succ, double* objval);

/**
 Produce the temperature list.
 
 @param inst instance of the struct "instance" for TSP problem.
 @param temperature_list temperature list.
 @param temperature_size temperature list size.
 */
void produce_temperature_list(instance* inst, double* temperature_list, int temperature_size);

//////////////////////////// utils ////////////////////////////

/**
 Return the max temperature and the corresponding array index.
 
 @param temp_list temperature list.
 @param size temperature list size.
 @param t_max value of max temperature.
 @param index array index of temperature list with max temperature.
 */
void t_max(double* temp_list, int size, double* t_max, int* index);

/**
 @brief Compute the distance between two nodes.
 
 @param i first node.
 @param j second node.
 @param inst instance of the struct "instance" for TSP problem.
 @return distance between two nodes.
 */
double dist(int i, int j, instance *inst);

/**
 Reverse the direction of a part of the tour.
 
 @param inst instance of the struct "instance" for TSP problem.
 @param start first node.
 @param end last node.
 @param succ TSP tour described with successors.
 */
void reverse_segment(instance* inst, int start, int end, int* succ);

#endif /* simulated_annealing_h */
