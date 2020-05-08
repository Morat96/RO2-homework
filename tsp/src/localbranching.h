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

void localbranching(instance *inst, CPXENVptr env, CPXLPptr lp);

// utils
double second(void);
void print_error(const char *err);
int my_generic_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *user);

#endif /* localbranching_h */
