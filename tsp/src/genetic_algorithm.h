//
//  genetic_algorithm.h
//  cplex
//
//  Created by Matteo Moratello on 04/06/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//

#ifndef genetic_algorithm_h
#define genetic_algorithm_h

#include <stdio.h>
#include "tsp.h"

void genetic_algorithm(instance* inst, int generations, int population_size, int crossover_size);
void mutation(instance* inst, int* child);
void get_best_parents(instance* inst, int* parents, double* fitness, int parents_size, int pop_size);
void update_fitness(instance* inst, int** population, int old_size, int* old_pop, double* fitness);
void update_population(instance* inst, double* fitness, int size_pop, int* old_pop, int num_old_pop);
void get_random_parents(instance* inst, int* parents, int parents_size, int pop_size);
void shortcut(instance* inst, int* new_child);
void repair_child(instance* inst, int* child, int size_child, int* nic, int size_nic);
void crossover(instance* inst, int* parents, int size, int** population, int** childs, int pop_size);
void get_random_population(instance* inst, int** population, int size, double* fitness);
void compute_mean_best_fitness(double* fitness, int size, double* mean, double* best, int* best_index);
void get_random_population(instance* inst, int** population, int size, double* fitness);

// utils
int remove_index(int* index, int from, int to);
double dist(int i, int j, instance *inst);
int three_min(double *array, int arr_size, int *ind);
int min(double *array, int arr_size);

#endif /* genetic_algorithm_h */
