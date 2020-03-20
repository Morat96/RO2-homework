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

#define VERBOSE 50 // printing level  (=10 only incumbent, =20 little output, =50-60 good, =70 verbose, >=100 cplex log)

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
    
} instance;


#endif /* tsp_h */
