//
//  vns.h
//  cplex
//
//  Created by Matteo Moratello on 04/06/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//

#ifndef vns_h
#define vns_h

#include <stdio.h>
#include "tsp.h"

void vns(instance* inst, int iter, int k);

// utils
double dist(int i, int j, instance *inst);
void twOpt(instance* inst, int *succ, double* objval);
int remove_index(int* index, int from, int to);

#endif /* vns_h */
