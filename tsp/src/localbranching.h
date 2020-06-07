//
//  localbranching.h
//  cplex
//
//  Created by Matteo Moratello on 27/04/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//

#ifndef localbranching_h
#define localbranching_h

#include <stdio.h>
#include "tsp.h"

/**
 Local branching algorithm.
 
 @param inst instance of the struct "instance" for TSP problem.
 @param env CPLEX environment.
 @param lp CPLEX LP.
 */
void localbranching(instance *inst, CPXENVptr env, CPXLPptr lp);

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
 Generic callback implementation.
 
 @param context Pointer to a callback context.
 @param contextid context identification.
 @param user Pointer to private user data.
 @return The routine returns 0 (zero) if successful and nonzero if an error occurs.
 */
int my_generic_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *user);

#endif /* localbranching_h */
