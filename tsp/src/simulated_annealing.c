//
//  simulated_annealing.c
//  cplex
//
//  Created by Matteo Moratello on 04/06/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//

#include "simulated_annealing.h"

// *********************************** Simulated Annealing ********************************* //
//
// List-Based Simulated Annealing Algorithm for Traveling Salesman Problem
// Shi-hua Zhan, Juan Lin, Ze-jun Zhang and Yi-wen Zhong
// https://www.hindawi.com/journals/cin/2016/1712630/
//
// ***************************************************************************************** //

/**
 Simulated annealing algorithm for TSP.

 @param inst instance of the struct "instance" for TSP problem.
 @param iter number of iterations.
 @param size temperature list size.
 */
void simulated_annealing(instance* inst, int iter, int size) {
    
    // produce the temperature list
    double* temp_list = (double*) calloc(size, sizeof(double));
    produce_temperature_list(inst, temp_list, size);
    double temp_max = 0.0;
    int index = 0;
    
    // find a first solution with NN
    double* xstar = (double*) calloc(((inst -> nnodes * (inst -> nnodes - 1))/2), sizeof(double));
    
    NearNeigh(inst, xstar);
    
    int *succ = (int *) calloc(inst -> nnodes, sizeof(int));
    int *comp = (int *) calloc(inst -> nnodes, sizeof(int));
    int ncomp = 0;
    
    build_sol(xstar, inst, succ, comp, &ncomp);
    
    // save the current solution as best solution
    int* best_sol = (int*) calloc(inst -> nnodes, sizeof(int));
    double best_objval = 0.0;
    for (int i = 0; i < inst -> nnodes; i++) {
        best_sol[i] = succ[i];
        best_objval += dist(i, best_sol[i], inst);
    }
    
    double objval = 0.0;
    double prob = 0.0;
    
    // outer loop
    for (int k = 0; k < iter; k++) {
        // keep the max temperature of the list
        t_max(temp_list, size, &temp_max, &index);
        double temp = 0;
        int c = 0;
        // inner loop
        for (int m = 0; m < inst -> nnodes; m++) {
            for (int i = 0; i < inst -> nnodes; i++) succ[i] = best_sol[i];
            // random neighbour solution
            neigh_random_swap(inst, succ, &objval);
            // if it is a better solution, set as best solution
            if (objval < best_objval) {
                for (int i = 0; i < inst -> nnodes; i++) best_sol[i] = succ[i];
                best_objval = objval;
            } else {
                // compute the probability
                prob = exp(-((objval - best_objval)/temp_max));
                double r = (double) rand()/RAND_MAX;
                if (r < prob) {
                    temp = temp - ((objval - best_objval) / log(r));
                    c++;
                    for (int i = 0; i < inst -> nnodes; i++) best_sol[i] = succ[i];
                    best_objval = objval;
                }
            }
        }
        if (c != 0) {
            temp_list[index] = temp / c;
        }
    }
    
    printf("\nBest objective function value: %lf\n", best_objval);
    print_solution(inst, succ);
    free(xstar);
    free(comp);
    free(succ);
    free(best_sol);
    free(temp_list);
}

/**
 Produce the temperature list.

 @param inst instance of the struct "instance" for TSP problem.
 @param temperature_list temperature list.
 @param temperature_size temperature list size.
 */
void produce_temperature_list(instance* inst, double* temperature_list, int temperature_size) {
    
    double* xstar = (double*) calloc(((inst -> nnodes * (inst -> nnodes - 1))/2), sizeof(double));
    
    // Phase 1
    // find a solution with NN and find a local optima using 2-OPT
    NearNeigh(inst, xstar);
    //random_solution(inst, xstar);
    
    int *succ = (int *) calloc(inst -> nnodes, sizeof(int));
    int *comp = (int *) calloc(inst -> nnodes, sizeof(int));
    int ncomp = 0;
    
    build_sol(xstar, inst, succ, comp, &ncomp);
    
    double objval = 0.0;
    double best_objval = 0.0;
    int* best_sol = (int*) calloc(inst -> nnodes, sizeof(int));
    for (int i = 0; i < inst -> nnodes; i++) {
        best_sol[i] = succ[i];
        best_objval += dist(i, best_sol[i], inst);
        
    }
    
    double p0 = 0.7;
    int i = 0;
    
    while (i < temperature_size) {
        for (int i = 0; i < inst -> nnodes; i++) succ[i] = best_sol[i];
        neigh_random_swap(inst, succ, &objval);
        if (objval < best_objval) {
            for (int i = 0; i < inst -> nnodes; i++) best_sol[i] = succ[i];
            best_objval = objval;
        }
        double temperature = - (fabs(objval - best_objval)/log(p0));
        temperature_list[i] = temperature;
        i++;
    }
    
    free(xstar);
    free(succ);
    free(comp);
    free(best_sol);
}

/**
 Swap two random edges of TSP tour.

 @param inst instance of the struct "instance" for TSP problem.
 @param succ TPS solution with successors format.
 @param objval TSP objective function value.
 */
void neigh_random_swap(instance* inst, int* succ, double* objval) {
    
    (*objval) = 0.0;
    
    // choose two random vertices
    int nFirst = rand() % inst -> nnodes;
    int nSecond = nFirst;
    
    while (nFirst == nSecond) {
        nSecond = rand() % inst -> nnodes;
    }
    
    for (int i = 0; i < inst -> nnodes; i++) {
        if (succ[i] == nFirst) break;
        if (succ[i] == nSecond) {
            int temp = nSecond;
            nSecond = nFirst;
            nFirst = temp;
            break;
        }
    }
    
    int succ_f = succ[nFirst];
    int succ_s = succ[nSecond];
    if (succ_f == nSecond) {
        nSecond = succ[nSecond];
        succ_s = succ[nSecond];
    }
    
    reverse_segment(inst, succ_f, nSecond, succ);
    succ[nFirst] = nSecond;
    succ[succ_f] = succ_s;
    
    for (int i = 0; i < inst -> nnodes; i++) (*objval) += dist(i, succ[i], inst);
}

/**
 Return the max temperature and the corresponding array index.

 @param temp_list temperature list.
 @param size temperature list size.
 @param t_max value of max temperature.
 @param index array index of temperature list with max temperature.
 */
void t_max(double* temp_list, int size, double* t_max, int* index) {
    (*t_max) = 0.0;
    for (int i = 0; i < size; i++) {
        if (temp_list[i] > (*t_max)) {
            (*t_max) = temp_list[i];
            (*index) = i;
        }
    }
}
