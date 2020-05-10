//
//  tsp.h
//  cplex
//
//  Created by Matteo Moratello on 13/03/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//

#ifndef tsp_h
#define tsp_h

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <cplex.h>
#include <pthread.h>
#include <time.h>

#define VERBOSE 1 // printing level  (= 10 only incumbent, = 20 little output, = 50-60 good, = 70 verbose, >= 100 cplex log)
#define debug 1   // debug checkings

// Data structures
typedef struct {
    
    //input data
    int nnodes;                              // number of nodes
    double *xcoord;                          // x-axis coords
    double *ycoord;                          // y-axis coords
    
    //parameters
    int model_type;                          // model type [0,7]
    double timelimit;                        // overall time limit, in sec.s
    char input_file[1000];                   // input file
    int randomseed;                          // random seed
    int loop;                                // resolve instance using Loop method {0,1}
    int callback;                            // resolve instance using CPLEX callbacks {0,1}
    int hardfixing;                          // resolve instance using Hard Fixing method {0,1}
    int localbranching;                      // resolve instance using Local Branching method {0,1}
    
    // distance type
    int euc_2d;                              // Euclidean 2D distance as metric distance
    int att;                                 // Att Euclidean distance as metric distance
    int geo;                                 // Geographic distance as metric distance
    
    // model;
    int ncols;                               // number of cols
    int ystart;                              // number of x(*,*) variables
    int ustart;                              // number of y(*,*) variables
    int zstart;                              // number of u(*,*) variables
    
    // heuristic solutions (one for each thread)
    int flag[4];                             // if the solution at thread "index" is updated {0,1}
    int **sol_thread;                        // heuristic solutions
    
} instance;

// struct for Concorde mincut
typedef struct {
    
    instance *inst;                          // TSP instance
    CPXCENVptr env;                          // CPLEX environment
    void *cbdata;                            // callback data
    int wherefrom;                           // callback info
    int *useraction_p;                       // callback user data
    
} input;

// **************** COMPACT MODELS ************** //
// STSP
void build_model_0(instance *inst, CPXENVptr env, CPXLPptr lp);
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
// ********************************************** //

// **************** LOOP METHOD **************** //
void add_constraints(instance *inst, CPXENVptr env, CPXLPptr lp, int *succ, int *comp, int ncomp, int n);
void loop_method(instance *inst, CPXENVptr env, CPXLPptr lp, double t1);
// ********************************************** //

// **************** LAZY CALLBACK *************** //
// lazy callback: integer LP
int CPXPUBLIC mylazycallback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p);
// add SEC in integer solutions
int myseparation(instance *inst, double *xstar, CPXCENVptr env, void *cbdata, int wherefrom);
// user callback: relaxation
int CPXPUBLIC UserCutCallback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p);
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

// **************** HARD FIXING **************** //
void hardfixing(instance *inst, CPXENVptr env, CPXLPptr lp);
// ********************************************* //

// ************** LOCAL BRANCHING ************** //
void localbranching(instance *inst, CPXENVptr env, CPXLPptr lp);
// ********************************************* //

// ***************** HEURISTICS **************** //
void NearNeigh(instance *inst, double *xstar);
void grasp(instance *inst, double *xstar);
void insertion(instance *inst, double *xstar);
void insertion_ch(instance *inst, double *xstar);
// ********************************************* //

double second(void);
void print_error(const char *err);
void print_solution(instance *inst, int *succ);
void print_solution_light(instance *inst, int *succ);
void build_sol(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp);
void build_compact_sol(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp);

#endif /* tsp_h */
