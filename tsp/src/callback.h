//
//  callback.h
//  cplex
//
//  Created by Matteo Moratello on 25/04/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//

#ifndef callback_h
#define callback_h

#include <stdio.h>
#include "tsp.h"
#include <cut.h>

// print error
void print_error(const char *err);
// position
int xpos(int i, int j, instance *inst);
// distance
double dist(int i, int j, instance *inst);
// build solution
void build_sol(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp);

// **************** LAZY CALLBACK *************** //
// lazy callback: integer LP
int CPXPUBLIC mylazycallback(CPXCENVptr env,
                             void *cbdata,
                             int wherefrom,
                             void *cbhandle,
                             int *useraction_p);
// add SEC in integer solutions
int myseparation(instance *inst, double *xstar, CPXCENVptr env, void *cbdata, int wherefrom);
// user callback: relaxation
int CPXPUBLIC UserCutCallback(CPXCENVptr env,
                              void *cbdata,
                              int wherefrom,
                              void *cbhandle,
                              int *useraction_p);
// add SEC in continuous relaxation solutions
int doit_fn_concorde(double cutval , int cutcount , int *cut , void *inParam);
// add a tsp solution from CPLEX's integer solution
int CPXPUBLIC myheuristic (CPXCENVptr env,
                           void       *cbdata,
                           int        wherefrom,
                           void       *cbhandle,
                           double     *objval_p,
                           double     *x,
                           int        *checkfeas_p,
                           int        *useraction_p);
// ********************************************** //

// ************** GENERIC CALLBACK ************** //
// either integer and relaxation callback
int my_generic_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *user);
// add SEC in integer solutions
int my_separation(instance *inst, double *xstar, CPXCALLBACKCONTEXTptr context);
// add SEC in continuous relaxation solutions
int doit_fn_concorde_gen(double cutval, int cutcount, int *cut , void *inParam);
// ********************************************** //

// ************ MERGE COMP ALGORITHM ************ //
// build a cycle from the distinct components found by CPLEX
void complete_cycle(instance *inst, int *succ, int *comp, int *ncomp);
// ********************************************** //

#endif /* callback_h */
