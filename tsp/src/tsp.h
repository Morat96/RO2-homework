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

#define VERBOSE 1 // printing level  (=10 only incumbent, =20 little output, =50-60 good, =70 verbose, >=100 cplex log)
#define debug 0

//data structures
typedef struct {
    
    //input data
    int nnodes;
    double *demand;
    double *xcoord;
    double *ycoord;
    
    //parameters
    int model_type;
    double timelimit;                        // overall time limit, in sec.s
    char input_file[1000];                   // input file
    int randomseed;
    int loop;
    int callback;
    int hardfixing;
    int localbranching;
    
    // distance type
    int euc_2d;
    int att;
    int geo;
    
    // model;
    int ncols;
    int ystart;
    int ustart;
    int zstart;
    
} instance;

// struct for concorde mincut
typedef struct {
    
    instance *inst;
    CPXCENVptr env;
    void *cbdata;
    int wherefrom;
    int *useraction_p;
    
} input;

#endif /* tsp_h */
