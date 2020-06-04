//
//  multi_start.h
//  cplex
//
//  Created by Matteo Moratello on 04/06/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//

#ifndef multi_start_h
#define multi_start_h

#include <stdio.h>
#include "tsp.h"

void multi_start(instance* inst, int iter);

// utils
double dist(int i, int j, instance *inst);
void random_solution_ms(instance* inst, int* succ, double* objval);
void twOpt(instance* inst, int *succ, double* objval);

#endif /* multi_start_h */
