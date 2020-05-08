//
//  hardfixing.h
//  cplex
//
//  Created by Matteo Moratello on 26/04/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//

#ifndef hardfixing_h
#define hardfixing_h

#include <stdio.h>
#include "tsp.h"

double second(void);
void print_error(const char *err);
int my_generic_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *user);

void hardfixing(instance *inst, CPXENVptr env, CPXLPptr lp);

#endif /* hardfixing_h */
