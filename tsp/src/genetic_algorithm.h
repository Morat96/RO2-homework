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

/**
 Genetic algorithm for TSP.
 
 @param inst instance of the struct "instance" for TSP problem.
 @param generations number of generations.
 @param population_size size of population.
 @param crossover_size number of new individual to generate at each generation.
 */
void genetic_algorithm(instance* inst, int generations, int population_size, int crossover_size);

/**
 Compute the first random population, that is, array of TSP feasible tours.
 
 @param inst instance of the struct "instance" for TSP problem.
 @param population matrix of size (#individuals) * (# nodes) containing the population.
 @param size population size, that is number of individuals.
 @param fitness objective function values of the populations, ie tour lenghts.
 */
void get_random_population(instance* inst, int** population, int size, double* fitness);

/**
 Get parents for generate new generation of individuals.
 With this function we select parents that have a low objective function with higher probability, that is best TSP solutions.
 
 @param inst instance of the struct "instance" for TSP problem.
 @param parents indices of population array of selected parents.
 @param fitness array with objective function values of the population.
 @param parents_size number of parents to select.
 @param pop_size size of the population.
 */
void get_best_parents(instance* inst, int* parents, double* fitness, int parents_size, int pop_size);

/**
 Get parents for generate new generation of individuals.
 With this function we select parents at random.
 
 @param inst instance of the struct "instance" for TSP problem.
 @param parents indices of population array of selected parents.
 @param parents_size number of parents to select.
 @param pop_size size of the population.
 */
void get_random_parents(instance* inst, int* parents, int parents_size, int pop_size);

/**
 Compute new individuals from parents using crossover algorithm.
 
 @param inst instance of the struct "instance" for TSP problem.
 @param parents parents array.
 @param size number of new childs of next generation.
 @param population population matrix.
 @param childs matrix of new childs.
 @param pop_size population size.
 */
void crossover(instance* inst, int* parents, int size, int** population, int** childs, int pop_size);

/**
 Compute the new individual combining two parents.
 
 @param inst instance of the struct "instance" for TSP problem.
 @param new_child new individual, i.e. new feasible TSP tour.
 */
void shortcut(instance* inst, int* new_child);

/**
 Function that "repair" a child using insertion algorithm.
 
 @param inst instance of the struct "instance" for TSP problem.
 @param child first part of the new individual.
 @param size_child size of the first part of the individual.
 @param nic nodes not in the solution.
 @param size_nic number of nodes not in the solution.
 */
void repair_child(instance* inst, int* child, int size_child, int* nic, int size_nic);

/**
 Apply a mutation to an individual.
 
 @param inst instance of the struct "instance" for TSP problem.
 @param child individual, i.e. TSP tour.
 */
void mutation(instance* inst, int* child);

/**
 Select individuals to replace in the next generation.
 @brief This function select the individuals based on their objective function value, that is, lower is their objective function
 higher is the probability to been choose as replace.
 
 @param inst instance of the struct "instance" for TSP problem.
 @param fitness objective function values of the population.
 @param size_pop population size.
 @param old_pop indices of individuals to replace.
 @param num_old_pop number of individuals to replace.
 */
void update_population(instance* inst, double* fitness, int size_pop, int* old_pop, int num_old_pop);

/**
 Update the fitness array, that is, lenghts of population TSP tours.
 
 @param inst instance of the struct "instance" for TSP problem.
 @param population population matrix.
 @param old_size number of new childs.
 @param old_pop indices of population of old inviduals.
 @param fitness objective function values of population, i.e. TSP tour lenghts.
 */
void update_fitness(instance* inst, int** population, int old_size, int* old_pop, double* fitness);

/**
 Compute new mean and best objective function value.
 
 @param fitness instance of the struct "instance" for TSP problem.
 @param size number of individuals in the population.
 @param mean mean of the objective function values.
 @param best best objective function values.
 @param best_index index of the best objective function value.
 */
void compute_mean_best_fitness(double* fitness, int size, double* mean, double* best, int* best_index);

//////////////////////////// utils ////////////////////////////

/**
 Remove element from an array by sliding it.
 
 @param index array of int.
 @param from start index.
 @param to size of array.
 @return if removal was successful.
 */
int remove_index(int* index, int from, int to);

/**
 @brief Compute the distance between two nodes.
 
 @param i first node.
 @param j second node.
 @param inst instance of the struct "instance" for TSP problem.
 @return distance between two nodes.
 */
double dist(int i, int j, instance *inst);

/**
 Compute the three smaller values of an array and store the indexes in the "ind" array.
 
 @param array array of double.
 @param arr_size array size.
 @param ind array of int containing the three smallest values of array.
 @return return 1 if the array size is smaller than three, otherwise return 0.
 */
int three_min(double *array, int arr_size, int *ind);

/**
 Compute minimum of an array and return the corresponding index.
 
 @param array array of double.
 @param arr_size array size.
 @return array index of minimum value of the array.
 */
int min(double *array, int arr_size);

#endif /* genetic_algorithm_h */
