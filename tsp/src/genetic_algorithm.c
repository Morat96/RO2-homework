//
//  genetic_algorithm.c
//  cplex
//
//  Created by Matteo Moratello on 04/06/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//

#include "genetic_algorithm.h"

// ************************************ Genetic Algorithm ********************************** //
//
// Description: TO DO.
//
// ***************************************************************************************** //

/**
 Genetic algorithm for TSP.

 @param inst instance of the struct "instance" for TSP problem.
 @param generations number of generations.
 @param population_size size of population.
 @param crossover_size number of new individual to generate at each generation.
 */
void genetic_algorithm(instance* inst, int generations, int population_size, int crossover_size) {
    
    int** population = (int**) calloc(population_size, sizeof(int*));
    double* fitness = (double*) calloc(population_size, sizeof(double));
    double avg_fitness = 0.0;
    double best_fitness = 0.0;
    for (int i = 0; i < population_size; i++) population[i] = calloc(inst -> nnodes, sizeof(int));
    int* parents = (int*) calloc(crossover_size * 2, sizeof(int));
    int** childs = (int**) calloc(crossover_size, sizeof(int*));
    for (int i = 0; i < crossover_size; i++) childs[i] = calloc(inst -> nnodes, sizeof(int));
    int* old_pop = (int*) calloc(crossover_size, sizeof(int));
    int index_best_objval = 0;
    
    // first population using Grasp
    get_random_population(inst, population, population_size, fitness);
    
    for (int iter = 0; iter < generations; iter ++) {
        
        double start_iter = second();
        
        get_best_parents(inst, parents, fitness, crossover_size * 2, population_size);
        //get_random_parents(inst, parents, crossover_size * 2, population_size);
        
        crossover(inst, parents, crossover_size, population, childs, population_size);
        
        // update population
        update_population(inst, fitness, population_size, old_pop, crossover_size);
        for (int i = 0; i < crossover_size; i++) {
            int index = old_pop[i];
            for (int j = 0; j < inst -> nnodes; j++) {
                population[index][j] = childs[i][j];
            }
        }
        
        update_fitness(inst, population, crossover_size, old_pop, fitness);
        compute_mean_best_fitness(fitness, population_size, &avg_fitness, &best_fitness, &index_best_objval);
        
        double stop_iter = second();
        
        printf("Generation: %d\n", iter + 1);
        printf("Average objective function value: %lf\n", avg_fitness);
        printf("Best objective function value: %lf\n", best_fitness);
        printf("Time : %lf\n\n", stop_iter - start_iter);
    }
    
    printf("Best objective function value: %lf\n", fitness[index_best_objval]);
    
    int* succ = (int*) calloc(inst -> nnodes, sizeof(int));
    for (int i = 0; i < inst -> nnodes - 1; i++) {
        int curr = population[index_best_objval][i];
        succ[curr] = population[index_best_objval][i+1];
    }
    int curr = population[index_best_objval][inst -> nnodes - 1];
    succ[curr] = population[index_best_objval][0];
    
    print_solution_light(inst, succ);
    
    // free
    for (int i = 0; i < population_size; i++) free(population[i]);
    for (int i = 0; i < crossover_size; i++) free(childs[i]);
    free(childs);
    free(population);
    free(fitness);
    free(parents);
    free(old_pop);
    free(succ);
}

/**
 Compute the first random population, that is, array of TSP feasible tours.

 @param inst instance of the struct "instance" for TSP problem.
 @param population matrix of size (#individuals) * (# nodes) containing the population.
 @param size population size, that is number of individuals.
 @param fitness objective function values of the populations, ie tour lenghts.
 */
void get_random_population(instance* inst, int** population, int size, double* fitness) {
    
    // number of nodes of the instance
    int n = inst -> nnodes;
    
    int *indices = (int *) calloc(inst -> nnodes - 1, sizeof(int));
    double *distances = (double *) calloc(inst -> nnodes - 1, sizeof(double));
    int *ind = (int *) calloc(3, sizeof(int));
    
    for (int p = 0; p < size; p ++) {
        
        // counter of resolved nodes
        int cnt = 1;
        
        // first node [0, nnodes - 1]
        int start_node = 0;
        int n_ind= 0;
        // variable that store the current node (at the beginning the current node is the start node)
        int current_node = start_node;
        population[p][n_ind++] = start_node;
        
        int nind = 0;
        // indices of nodes except the current node
        for (int i = 0; i < n - cnt + 1; i++) if (i != current_node) indices[nind++] = i;
        
        // obj function value
        double objval = 0.0;
        
        // randomize seed
        //srand((unsigned int)time(NULL));
        
        // in each cycle find the minimum edge
        while (cnt < inst -> nnodes - 2) {
            
            // compute distances
            for (int i = 0; i < n - cnt; i++) distances[i] = dist(current_node, indices[i], inst);
            
            // find the three shorter edges from the current node
            int minim = 0;
            if (three_min(distances, n - cnt, ind)) print_error("Error in function three_min");
            
            /* choose the next node with weighted probability.
             The node corresponding to the smaller edge has 50% probability to be chosen.
             The other two nodes have 25% probability to be chosen. */
            
            // generate a random value from [0,1]
            float rv = ((float)rand()/(float)(RAND_MAX));
            
            // set the new edge based on probabilities
            if (rv <= 0.34) minim = ind[0];
            else if (rv <= 0.67) minim = ind[1];
            else minim = ind[2];
            
            // update objective function value
            objval += distances[minim];
            
            // update the solution
            population[p][n_ind++] = indices[minim];
            
            // update the current node
            current_node = indices[minim];
            
            // remove the previous node
            for (int c = minim; c < n - cnt; c++) {
                distances[c] = distances[c+1];
                indices[c] = indices[c+1];
            }
            
            // shrink the array of one element
            cnt ++;
        }
        
        // choose the last two edges with the greedy method
        // only two and one edge to choose from
        for (int i = 0; i < 2; i++) {
            
            // compute distances
            for (int i = 0; i < n - cnt; i++) distances[i] = dist(current_node, indices[i], inst);
            
            // find the shorter edge from the current node
            int minim = min(distances, n - cnt);
            
            // update objective function value
            objval += distances[minim];
            
            // update the solution
            population[p][n_ind++] = indices[minim];
            
            // update the current node
            current_node = indices[minim];
            
            // remove the previous node
            for (int c = minim; c < n - cnt; c++) {
                distances[c] = distances[c+1];
                indices[c] = indices[c+1];
            }
            
            // shrink the array of one element
            cnt ++;
        }
        
        objval += dist(current_node, start_node, inst);
        
        //printf("Objective function value: %lf\n", objval);
        fitness[p] = objval;
    }
    
    free(ind);
    free(distances);
    free(indices);
}

/**
 Get parents for generate new generation of individuals.
 With this function we select parents that have a low objective function with higher probability, that is best TSP solutions.

 @param inst instance of the struct "instance" for TSP problem.
 @param parents indices of population array of selected parents.
 @param fitness array with objective function values of the population.
 @param parents_size number of parents to select.
 @param pop_size size of the population.
 */
void get_best_parents(instance* inst, int* parents, double* fitness, int parents_size, int pop_size) {
    
    double* cost = (double*) calloc(pop_size, sizeof(double));
    int count = (fitness[0] == 0) ? 1  : (log10(fitness[0]) + 1);
    for (int i = 0; i < pop_size; i++) cost[i] = (1/fitness[i])*(pow(10, count));
    
    double sum = 0.0;
    for (int i = 0; i < pop_size; i++) sum += cost[i];
    
    int num = 0;
    
    while (num != parents_size) {
        float rv = (float) rand()/RAND_MAX;
        double current_sum = 0.0;
        int cnt = 0;
        while (current_sum < rv) {
            current_sum += cost[cnt++]/sum;
        }
        int not_present = 1;
        for (int i = 0; i < num; i++) {
            if ((cnt - 1) == parents[i]) not_present = 0;
        }
        if (not_present) parents[num++] = cnt - 1;
    }
    
    free(cost);
}

/**
 Get parents for generate new generation of individuals.
 With this function we select parents at random.
 
 @param inst instance of the struct "instance" for TSP problem.
 @param parents indices of population array of selected parents.
 @param parents_size number of parents to select.
 @param pop_size size of the population.
 */
void get_random_parents(instance* inst, int* parents, int parents_size, int pop_size) {
    
    int* indeces = (int*) calloc(pop_size, sizeof(int));
    for (int i = 0; i < pop_size; i++) indeces[i] = i;
    int cnt = 0;
    
    while (cnt < parents_size) {
        int rv = rand() % (pop_size - cnt);
        parents[cnt] = indeces[rv];
        remove_index(indeces, rv, pop_size - cnt);
        cnt ++;
    }
    
    free(indeces);
}

/**
 Compute new individuals from parents using crossover algorithm.

 @param inst instance of the struct "instance" for TSP problem.
 @param parents parents array.
 @param size number of new childs of next generation.
 @param population population matrix.
 @param childs matrix of new childs.
 @param pop_size population size.
 */
void crossover(instance* inst, int* parents, int size, int** population, int** childs, int pop_size) {
    
    int* new_child = (int*) calloc(inst -> nnodes, sizeof(int));
    
    for (int i = 0; i < size; i++) {
        // compute the threshold for crossover
        int threshold = rand() % inst -> nnodes;
        int parent1 = parents[i * 2];
        int parent2 = parents[i * 2 + 1];
        for (int j = 0; j < threshold; j++) new_child[j] = population[parent1][j];
        for (int j = threshold; j < inst -> nnodes; j++) new_child[j] = population[parent2][j];
        //printf("|%d|\n", threshold);
        shortcut(inst, new_child);
        float rv = (float) rand()/RAND_MAX;
        if (rv < 0.1) mutation(inst, new_child);
        for (int add = 0; add < inst -> nnodes; add++) childs[i][add] = new_child[add];
    }
    
    free(new_child);
}

/**
 Compute the new individual combining two parents.

 @param inst instance of the struct "instance" for TSP problem.
 @param new_child new individual, i.e. new feasible TSP tour.
 */
void shortcut(instance* inst, int* new_child) {
    
    int* rep = (int*) calloc(inst -> nnodes, sizeof(int));
    int* not_in_child = (int*) calloc(inst -> nnodes, sizeof(int));
    
    int size = 0;
    int cnt = 0;
    int cnt2 = 0;

    for (int i = 0; i < inst -> nnodes; i++) {
        int nic = 0;
        for (int j = 0; j < inst -> nnodes - cnt; j++) {
            if (new_child[j] == i) {
                for (int z = 0; z < size; z++) {
                    if (new_child[j] == rep[z]) {
                        for (int m = j; m < inst -> nnodes; m++) new_child[m] = new_child[m + 1];
                        cnt ++;
                    }
                }
                rep[size++] = i;
            }
            else nic++;
        }
        if (nic == inst -> nnodes - cnt) not_in_child[cnt2++] = i;
    }
    
    repair_child(inst, new_child, inst -> nnodes - cnt, not_in_child, cnt2);
    
    free(rep);
    free(not_in_child);
    
}

/**
 Function that "repair" a child using insertion algorithm.

 @param inst instance of the struct "instance" for TSP problem.
 @param child first part of the new individual.
 @param size_child size of the first part of the individual.
 @param nic nodes not in the solution.
 @param size_nic number of nodes not in the solution.
 */
void repair_child(instance* inst, int* child, int size_child, int* nic, int size_nic) {
    
    // number of nodes
    int n = inst -> nnodes;
    
    int start = 0;
    int n_nodes_sol = size_child;
    
    // End first part
    // Second part
    // insert nodes on the solution until all edges belong to C
    
    int cnt = size_child;
    double objval = 0.0;
    
    // cycle until all edges belong to C
    while (n_nodes_sol != n) {
        
        // variables for computing the best edge
        double best_cost = INT_MAX;
        int best_h = 0;
        int best_pos = 0;
        int best_h_index = 0;
        
        // for each edge in C
        for (int c = 0; c < n_nodes_sol - 1; c++) {
            // for each remaining node
            for (int i = 0; i < n - cnt; i++) {
                // compute the cost function
                // c_ah + c_hb - c_ab
                double cost_h = dist(child[c], nic[i], inst) + dist(nic[i], child[c + 1], inst) - dist(child[c], child[c + 1], inst);
                
                // updating the best cost value
                if (cost_h < best_cost) {
                    best_h = nic[i];
                    best_h_index = i;
                    best_pos = c;
                    best_cost = cost_h;
                }
            }
        }
        
        // last edge (the one that closes the cycle)
        for (int i = 0; i < n - cnt; i++) {
            
            double cost_h = dist(child[n_nodes_sol - 1], nic[i], inst) + dist(nic[i], child[start], inst) - dist(child[n_nodes_sol - 1], child[start], inst);
            
            // updating the best cost value
            if (cost_h < best_cost) {
                best_h = nic[i];
                best_h_index = i;
                best_pos = n_nodes_sol - 1;
                best_cost = cost_h;
            }
        }
        
        // update the solution with the best edge
        for (int c = n - 1; c > best_pos; c--) child[c] = child[c - 1];
        child[best_pos + 1] = best_h;
        
        // remove the previous node
        for (int c = best_h_index; c < n - cnt; c++) nic[c] = nic[c + 1];
        
        cnt ++;
        n_nodes_sol ++;
    }
    
    // build the xstar and the objective function value
    for (int i = 0; i < n - 1; i++) {
        objval += dist(child[i], child[i+1], inst);
    }
    
    objval += dist(child[n - 1], child[start], inst);
    
    //printf("Objective function value: %lf\n", objval);
    
}

/**
 Apply a mutation to an individual.

 @param inst instance of the struct "instance" for TSP problem.
 @param child individual, i.e. TSP tour.
 */
void mutation(instance* inst, int* child) {
    
    int first_index = 0;
    int second_index = 0;
    
    while (first_index != second_index) {
        first_index = rand() % inst -> nnodes;
        second_index = rand() % inst -> nnodes;
    }
    int temp = child[first_index];
    child[first_index] = child[second_index];
    child[second_index] = temp;
}

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
void update_population(instance* inst, double* fitness, int size_pop, int* old_pop, int num_old_pop) {
    
    // compute probabilities
    double* prob = (double*) calloc(size_pop, sizeof(double));
    double sum = 0.0;
    for (int i = 0; i < size_pop; i++) sum += fitness[i];
    for (int i = 0; i < size_pop; i++) {
        prob[i] = fitness[i]/sum;
    }
    int num = 0;
    
    while (num != num_old_pop) {
        float rv = (float) rand()/RAND_MAX;
        double current_sum = 0.0;
        int cnt = 0;
        while (current_sum < rv) {
            current_sum += prob[cnt++];
        }
        int not_present = 1;
        for (int i = 0; i < num; i++) {
            if ((cnt - 1) == old_pop[i]) not_present = 0;
        }
        if (not_present) old_pop[num++] = cnt - 1;
    }
    
    free(prob);
}

/**
 Update the fitness array, that is, lenghts of population TSP tours.

 @param inst instance of the struct "instance" for TSP problem.
 @param population population matrix.
 @param old_size number of new childs.
 @param old_pop indices of population of old inviduals.
 @param fitness objective function values of population, i.e. TSP tour lenghts.
 */
void update_fitness(instance* inst, int** population, int old_size, int* old_pop, double* fitness) {
    
    for (int i = 0; i < old_size; i++) {
        double objval = 0.0;
        int index_pop = old_pop[i];
        for (int j = 0; j < inst -> nnodes - 1; j++) objval += dist(population[index_pop][j], population[index_pop][j+1], inst);
        objval += dist(population[index_pop][inst -> nnodes - 1], population[index_pop][0], inst);
        fitness[index_pop] = objval;
    }
}

/**
 Compute new mean and best objective function value.

 @param fitness instance of the struct "instance" for TSP problem.
 @param size number of individuals in the population.
 @param mean mean of the objective function values.
 @param best best objective function values.
 @param best_index index of the best objective function value.
 */
void compute_mean_best_fitness(double* fitness, int size, double* mean, double* best, int* best_index) {
    (*mean) = 0.0;
    (*best) = INT_MAX;
    for (int i = 0; i < size; i++) {
        (*mean) += fitness[i];
        if (fitness[i] < (*best)) {
            (*best) = fitness[i];
            (*best_index) = i;
        }
    }
    (*mean) /= size;
}





