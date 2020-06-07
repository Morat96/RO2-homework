//
//  compact_models.h
//  cplex
//
//  Created by Matteo Moratello on 13/04/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//

#ifndef compact_models_h
#define compact_models_h

#include <stdio.h>
#include "tsp.h"

// COMPACT MODELS

// compact model: MTZ
void build_model_1(instance *inst, CPXENVptr env, CPXLPptr lp);
// compact model: FLOW1
void build_model_2(instance *inst, CPXENVptr env, CPXLPptr lp);
// compact model: FLOW2
void build_model_3(instance *inst, CPXENVptr env, CPXLPptr lp);
// compact model: FLOW3
void build_model_4(instance *inst, CPXENVptr env, CPXLPptr lp);
// compact model: T1
void build_model_5(instance *inst, CPXENVptr env, CPXLPptr lp);
// compact model: T2
void build_model_6(instance *inst, CPXENVptr env, CPXLPptr lp);
// compact model: T3
void build_model_7(instance *inst, CPXENVptr env, CPXLPptr lp);

/**
 Build succ() and comp() wrt xstar().
 
 @param xstar solution in CPLEX format.
 @param inst instance of the struct "instance" for TSP problem.
 @param succ solution as successors.
 @param comp number of the component associated to each node in the solution.
 @param ncomp number of components in the solution.
 */
void build_compact_sol(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp);

// Position of variables of the CPLEX solutions.

/**
 For compact models solutions.
 
 @param i first node.
 @param j second node.
 @param inst instance of the struct "instance" for TSP problem.
 @return position on the CPLEX solution.
 */
int xpos_compact(int i, int j, instance *inst) { return i * inst->nnodes + j; }

/**
 For y(*,*) variables.

 @param i first node.
 @param j second node.
 @param inst instance of the struct "instance" for TSP problem.
 @return position on the CPLEX solution.
 */
int ypos(int i, int j, instance *inst) { return inst-> ystart + i * inst->nnodes + j; }

/**
 For y(*,*,*) variables.

 @param i first solution index.
 @param j second solution index.
 @param k third solution index.
 @param inst instance of the struct "instance" for TSP problem.
 @return position on the CPLEX solution.
 */
int y2pos(int i, int j, int k, instance *inst) { return inst-> ystart + i * inst->nnodes + j + (k-1) * (inst -> nnodes) * (inst -> nnodes); }

/**
 For y(*,*,*) variables.
 
 @param i first solution index.
 @param j second solution index.
 @param k third solution index.
 @param inst instance of the struct "instance" for TSP problem.
 @return position on the CPLEX solution.
 */
int y3pos(int i, int j, int k, instance *inst) { return inst-> ystart + i * inst->nnodes + j + (k) * (inst -> nnodes) * (inst -> nnodes); }

/**
 For z(*,*) variables.
 
 @param i first solution index.
 @param j second solution index.
 @param inst instance of the struct "instance" for TSP problem.
 @return position on the CPLEX solution.
 */
int zpos(int i, int j, instance *inst) { return inst-> zstart + i * inst->nnodes + j; }

/**
 For u(*) variables.

 @param i first solution index.
 @param inst instance of the struct "instance" for TSP problem.
 @return position on the CPLEX solution.
 */
int upos(int i, instance *inst) { return inst -> ustart + i; }

//////////////////////////// utils ////////////////////////////

/**
 Print an error and exit from the program.
 
 @param err error to show.
 */
void print_error(const char *err);

/**
 @brief Compute the distance between two nodes.
 
 @param i first node.
 @param j second node.
 @param inst instance of the struct "instance" for TSP problem.
 @return distance between two nodes.
 */
double dist(int i, int j, instance *inst);

#endif /* compact_models_h */
