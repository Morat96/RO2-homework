//
//  tabu_search.h
//  cplex
//
//  Created by Matteo Moratello on 04/06/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//

#ifndef tabu_search_h
#define tabu_search_h

#include <stdio.h>
#include "tsp.h"

typedef struct {
    int vert1, vert2;
} tabu_vert;

void tabu_search(instance* inst, int iter, int list_size);
void findBestEdgesSwap(instance* inst, int* succ, double* objval, tabu_vert* tabu_list, int* tabu_size, int list_size, int check);

// utils
double dist(int i, int j, instance *inst);
void twOpt(instance* inst, int *succ, double* objval);

#endif /* tabu_search_h */
