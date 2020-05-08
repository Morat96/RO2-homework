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

typedef struct tPoint {
    double x, y;
} Point;

int convexHull(Point points[], int n, Point* ch, int* size);

void NearNeigh(instance *inst, double *xstar);
void grasp(instance *inst, double *xstar);
void insertion(instance *inst, double *xstar);
void insertion_ch(instance *inst, double *xstar);

// utils
double dist(int i, int j, instance *inst);
int xpos(int i, int j, instance *inst);
void build_sol(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp);
void print_solution(instance *inst, int *succ);
void print_solution_light(instance *inst, int *succ);
void print_error(const char *err);
double second(void);

#endif /* heuristics_h */
