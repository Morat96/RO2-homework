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
    
} instance;

// struct for Concorde mincut
typedef struct {
    
    instance *inst;                          // TSP instance
    CPXCENVptr env;                          // CPLEX environment
    void *cbdata;                            // callback data
    int wherefrom;                           // callback info
    int *useraction_p;                       // callback user data
    
} input;

#endif /* tsp_h */
