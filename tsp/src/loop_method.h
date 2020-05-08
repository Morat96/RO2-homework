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

// loop method
void add_constraints(instance *inst, CPXENVptr env, CPXLPptr lp, int *succ, int *comp, int ncomp, int n);
void loop_method(instance *inst, CPXENVptr env, CPXLPptr lp, double t1);

// utils
double second(void);
void print_error(const char *err);
void print_solution(instance *inst, int *succ);
void build_sol(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp);
int xpos(int i, int j, instance *inst);

#endif /* loop_method_h */
