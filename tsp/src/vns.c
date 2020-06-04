//
//  vns.c
//  cplex
//
//  Created by Matteo Moratello on 04/06/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//

#include "vns.h"

// ************************** Variable Neighborhood Search ************************* //
//
// Description:
// - Build a first solution using NN algorithm and apply to it the 2-OPT algorithm.
// - Compute and apply a random k-move to the solution.
// - Reach the local minimum of the tour.
// - Repeat "iter" times and save the best solution found.
//
// Parameters:
// instance* inst: instance of the problem.
// int iter: number of iterations.
// int k: value for the random k-move.
//
// ********************************************************************************* //

void vns(instance* inst, int iter, int k) {
    
    if (k % 2 == 0) k += 1;
    
    double best_obj = 0.0;
    double current_obj = 0.0;
    // incumbent
    int *succ_best = (int *) calloc(inst -> nnodes, sizeof(int));
    
    // first solution with nearest neighborhood algorithm
    double* xstar = (double*) calloc(((inst -> nnodes * (inst -> nnodes - 1))/2), sizeof(double));
    NearNeigh(inst, xstar);
    
    int *succ = (int *) calloc(inst -> nnodes, sizeof(int));
    int *comp = (int *) calloc(inst -> nnodes, sizeof(int));
    int ncomp = 0;
    
    build_sol(xstar, inst, succ, comp, &ncomp);
    
    // search in a neighborhood with 2-OPT algorithm
    twOpt(inst, succ, &current_obj);
    best_obj = current_obj;
    // best solution found
    for (int i = 0; i < inst -> nnodes; i++) succ_best[i] = succ[i];
    
    int* index = (int*) calloc(inst -> nnodes, sizeof(int));
    int* rand_index = (int*) calloc(k, sizeof(int));
    int* ord_index = (int*) calloc(k, sizeof(int));
    int rv = 0;
    
    // algorithm iterations
    for (int i = 0; i < iter; i++) {
        
        // current object function value
        current_obj = 0.0;
        
        // k-OPT random move
        
        // select k distinct random indexes
        for (int i = 0; i < inst -> nnodes; i++) index[i] = i;
        for (int i = 0; i < k; i++) rand_index[i] = -1;
        rv = rand() % inst -> nnodes;
        rand_index[0] = index[rv];
        // remove the previous node
        remove_index(index, rv, inst -> nnodes);
        
        int cnt = 1;
        
        // k distinct random indices
        while (cnt != k) {
            
            int rv = rand() % (inst -> nnodes - cnt);
            rand_index[cnt] = index[rv];
            remove_index(index, rv, inst -> nnodes - cnt);
            cnt ++;
        }
        
        int cnt2 = 0;
        int n = 0;
        
        // select indices w.r.t. tour order
        for (int i = 0; n < inst -> nnodes; n++, i = succ[i]) {
            for (int j = 0; j < k - cnt2; j++) {
                if (i == rand_index[j]) {
                    ord_index[cnt2] = rand_index[j];
                    remove_index(rand_index, j, k - cnt2);
                    cnt2 ++;
                }
            }
        }
        
        // k-OPT random move
        int first = succ[ord_index[0]];
        for (int i = 0; i < k - 1; i++) {
            succ[ord_index[i]] = succ[ord_index[i + 1]];
        }
        succ[ord_index[k - 1]] = first;
        
        // intensification with 2-OPT
        twOpt(inst, succ, &current_obj);
        
        // check if it is a better solution than the best found so far
        if (current_obj < best_obj) {
            best_obj = current_obj;
            // update incumbent
            for (int i = 0; i < inst -> nnodes; i++) succ_best[i] = succ[i];
        }
    }
    
    printf("Best objective function value: %lf\n", best_obj);
    //print_solution_light(inst, succ_best);
    
    free(index);
    free(rand_index);
    free(ord_index);
    free(succ);
    free(succ_best);
    free(comp);
    free(xstar);
    
}
