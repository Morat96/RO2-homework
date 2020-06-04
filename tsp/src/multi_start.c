//
//  multi_start.c
//  cplex
//
//  Created by Matteo Moratello on 04/06/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//

#include "multi_start.h"

// *********************************** Multi Start ********************************* //
//
// Algorithm description:
// - Produce a random solution;
// - Reach the local minimum using 2-OPT algorithm;
// - Repeat "iter" times and save the best solution found.
//
// Parameters:
// instance* inst: instance of the problem.
// int iter: number of iterations.
//
// ********************************************************************************* //

void multi_start(instance* inst, int iter) {
    
    int* succ = (int*) calloc(inst -> nnodes, sizeof(int));
    int* best_succ = (int*) calloc(inst -> nnodes, sizeof(int));
    double objval = 0.0;
    double best_objval = INT_MAX;
    
    for (int i = 0; i < iter; i++) {
        random_solution_ms(inst, succ, &objval);
        twOpt(inst, succ, &objval);
        if (objval < best_objval) {
            best_objval = objval;
            for (int i = 0; i < inst -> nnodes; i++) best_succ[i] = succ[i];
        }
    }
    
    printf("Best objective function value: %lf\n", best_objval);
    print_solution_light(inst, best_succ);
    
    free(succ);
    free(best_succ);
}

// build a TSP random solution
void random_solution_ms(instance* inst, int* succ, double* objval) {
    
    int* sol = (int*) calloc(inst -> nnodes, sizeof(int));
    int* index = (int*) calloc(inst -> nnodes, sizeof(int));
    
    for (int i = 0; i < inst -> nnodes; i++) index[i] = i;
    int cnt = inst -> nnodes;
    int cnt2 = 0;
    
    while (cnt > 0) {
        
        int rv = rand() % cnt;
        
        sol[cnt2] = index[rv];
        
        // remove the previous node
        for (int c = rv; c < inst -> nnodes - cnt2; c++) {
            index[c] = index[c+1];
        }
        
        cnt2++;
        cnt --;
    }
    
    (*objval) = 0.0;
    
    for (int i = 0; i < inst -> nnodes - 1; i++) {
        int curr = sol[i];
        succ[curr] = sol[i + 1];
        (*objval) += dist(sol[i], sol[i + 1], inst);
    }
    int curr = sol[inst -> nnodes - 1];
    succ[curr] = sol[0];
    (*objval) += dist(sol[inst -> nnodes - 1], sol[0], inst);
    
    free(sol);
    free(index);
}
