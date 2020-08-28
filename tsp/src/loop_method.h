//
//  loop_method.h
//  cplex
//
//  Created by Matteo Moratello on 13/04/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//

#ifndef loop_method_h
#define loop_method_h

#include <stdio.h>
#include "tsp.h"

void add_constraints(instance *inst, CPXENVptr env, CPXLPptr lp, int *succ, int *comp, int ncomp, int n);

/**
 Loop method second version.
 
 @param inst instance of the struct "instance" for TSP problem.
 @param env CLEX environment.
 @param lp CLEX LP.
 */
void loop_method(instance *inst, CPXENVptr env, CPXLPptr lp, double t1);

void loop_method_vers1(instance *inst, CPXENVptr env, CPXLPptr lp);

//////////////////////////// utils ////////////////////////////

/**
 Compute the seconds passed from the program start.
 
 @return time in seconds.
 */
double second(void);

/**
 Print an error and exit from the program.
 
 @param err error to show.
 */
void print_error(const char *err);

/**
 Print the tour.
 
 @param inst instance of the struct "instance" for TSP problem.
 @param succ tour solution described with successors.
 */
void print_solution(instance *inst, int *succ);

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
 Return the index of CPLEX solution array from two subsequent nodes..
 
 @param i first node.
 @param j second node.
 @param inst instance of the struct "instance" for TSP problem.
 @return array position.
 */
int xpos(int i, int j, instance *inst);

#endif /* loop_method_h */
