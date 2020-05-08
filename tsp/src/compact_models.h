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

// BUILD COMPACT MODELS SOLUTION
void build_compact_sol(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp);

// utils
void print_error(const char *err);
void build_compact_sol(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp);

// variable position
int xpos_compact(int i, int j, instance *inst) { return i * inst->nnodes + j; }
int ypos(int i, int j, instance *inst) { return inst-> ystart + i * inst->nnodes + j; }
int y2pos(int i, int j, int k, instance *inst) { return inst-> ystart + i * inst->nnodes + j + (k-1) * (inst -> nnodes) * (inst -> nnodes); }
int y3pos(int i, int j, int k, instance *inst) { return inst-> ystart + i * inst->nnodes + j + (k) * (inst -> nnodes) * (inst -> nnodes); }
int zpos(int i, int j, instance *inst) { return inst-> zstart + i * inst->nnodes + j; }
int upos(int i, instance *inst) { return inst -> ustart + i; }

// distance between two nodes
double dist(int i, int j, instance *inst);

#endif /* compact_models_h */
