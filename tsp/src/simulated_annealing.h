//
//  simulated_annealing.h
//  cplex
//
//  Created by Matteo Moratello on 04/06/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//

#ifndef simulated_annealing_h
#define simulated_annealing_h

#include <stdio.h>
#include "tsp.h"

void simulated_annealing(instance* inst, int iter, int size);
void neigh_random_swap(instance* inst, int* succ, double* objval);
void produce_temperature_list(instance* inst, double* temperature_list, int temperature_size);

// utils
void t_max(double* temp_list, int size, double* t_max, int* index);
double dist(int i, int j, instance *inst);
void reverse_segment(instance* inst, int start, int end, int* succ);

#endif /* simulated_annealing_h */
