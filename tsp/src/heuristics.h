//
//  heuristics.h
//  cplex
//
//  Created by Matteo Moratello on 04/05/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//

#ifndef heuristics_h
#define heuristics_h

#include <stdio.h>
#include "tsp.h"

/**
 Point structure.
 */
typedef struct tPoint {
    double x, y;
} Point;

/**
 Compute the Convex Hull of a graph.

 @param points Array of Points (nodes of the graph).
 @param n number of nodes.
 @param ch Array of points of the Convex Hull.
 @param size number of points of the convex hull.
 @return 0 if successful, 1 otherwise.
 */
int convexHull(Point points[], int n, Point* ch, int* size);

/**
 Nearest Neighborhood (greedy).
 
 @param inst instance of the struct "instance" for TSP problem.
 @param xstar TSP solution using CPLEX format.
 */
void NearNeigh(instance *inst, double *xstar);

/**
 Grasp (randomization).
 
 @param inst instance of the struct "instance" for TSP problem.
 @param xstar TSP solution using CPLEX format.
 */
void grasp(instance *inst, double *xstar);

/**
 Insertion heuristic.
 
 @param inst instance of the struct "instance" for TSP problem.
 @param xstar TSP solution using CPLEX format.
 */
void insertion(instance *inst, double *xstar);

/**
 Insertion heuristic with Convex Hull.
 
 @param inst instance of the struct "instance" for TSP problem.
 @param xstar TSP solution using CPLEX format.
 */
void insertion_ch(instance *inst, double *xstar);

/**
 Refining algorithm: 2-OPT move.
 
 @param inst instance of the struct "instance" for TSP problem.
 @param succ TPS solution with successors format.
 @param objval TSP objective function value.
 */
void twOpt(instance* inst, int *succ, double* objval);
void twOptv2(instance* inst, double* xstar);

/**
 Refining algorithm: 3-OPT move.
 
 @param inst instance of the struct "instance" for TSP problem.
 @param xstar TSP solution using CPLEX format.
 */
void threeOpt(instance* inst, double* xstar);

/**
 build a TSP random solution
 
 @param inst instance of the struct "instance" for TSP problem.
 @param xstar TSP solution using CPLEX format.
 */
void random_solution(instance* inst, double* xstar);

//////////////////////////// utils ////////////////////////////

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

/**
 Reverse the direction of a part of the tour.
 
 @param inst instance of the struct "instance" for TSP problem.
 @param start first node.
 @param end last node.
 @param succ TSP tour described with successors.
 */
void reverse_segment(instance* inst, int start, int end, int* succ);

/**
 @brief Compute the distance between two nodes.
 
 @param i first node.
 @param j second node.
 @param inst instance of the struct "instance" for TSP problem.
 @return distance between two nodes.
 */
double dist(int i, int j, instance *inst);

/**
 Return the index of CPLEX solution array from two subsequent nodes..
 
 @param i first node.
 @param j second node.
 @param inst instance of the struct "instance" for TSP problem.
 @return array position.
 */
int xpos(int i, int j, instance *inst);

/**
 Build succ() and comp() wrt xstar().
 
 @param xstar CPLEX solution.
 @param inst instance of the struct "instance" for TSP problem.
 @param succ TSP solution as successors.
 @param comp component associated for each nodes.
 @param ncomp number of components in the solution.
 */
void build_sol(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp);

/**
 Print the tour.
 
 @param inst instance of the struct "instance" for TSP problem.
 @param succ tour solution described with successors.
 */
void print_solution(instance *inst, int *succ);

/**
 Print a light version of the tour.
 
 @param inst instance of the struct "instance" for TSP problem.
 @param succ tour solution described with successors.
 */
void print_solution_light(instance *inst, int *succ);

/**
 Print an error and exit from the program.

 @param err error to show.
 */
void print_error(const char *err);

/**
 Compute the seconds passed from the program start.

 @return time in seconds.
 */
double second(void);

/**
 Compute the k nearest nodes of each node in the graph.

 @param inst instance of the struct "instance" for TSP problem.
 @param distances matrix of nodes.
 */
void smallerKnodes(instance* inst, int** distances);

#endif /* heuristics_h */
