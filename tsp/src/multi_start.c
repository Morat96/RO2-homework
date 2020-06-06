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
// ********************************************************************************* //

/**
 Multi-start algorithm for TSP.

 @param inst instance of the struct "instance" for TSP problem.
 @param iter number of iterations.
 */
void multi_start(instance* inst, int iter) {
    
    int* succ = (int*) calloc(inst -> nnodes, sizeof(int));
    int* best_succ = (int*) calloc(inst -> nnodes, sizeof(int));
    double objval = 0.0;
    double best_objval = INT_MAX;
    
    for (int i = 0; i < iter; i++) {
        //random_solution_ms(inst, succ, &objval);
        grasp_solution_ms(inst, succ, &objval);
        twOpt(inst, succ, &objval);
        if (objval < best_objval) {
            best_objval = objval;
            for (int i = 0; i < inst -> nnodes; i++) best_succ[i] = succ[i];
        }
    }
    
    printf("Best objective function value: %lf\n", best_objval);
    //print_solution_light(inst, best_succ);
    
    free(succ);
    free(best_succ);
}

/**
 Build a random TSP feasible solution.

 @param inst instance of the struct "instance" for TSP problem.
 @param succ TSP tour represented with successors.
 @param objval TSP objective function value.
 */
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

/**
 Build a random TSP feasible solution using Grasp.

 @param inst inst instance of the struct "instance" for TSP problem.
 @param succ succ TSP tour represented with successors.
 @param objval objval TSP objective function value.
 */
void grasp_solution_ms(instance *inst, int* succ, double* objval) {
    
    (*objval) = 0.0;
    
    // number of nodes of the instance
    int n = inst -> nnodes;
    
    // counter of resolved nodes
    int cnt = 1;
    
    int *indices = (int *) calloc(inst -> nnodes - 1, sizeof(int));
    double *distances = (double *) calloc(inst -> nnodes - 1, sizeof(double));
    int *ind = (int *) calloc(3, sizeof(int));
    
    // first node [0, nnodes - 1]
    int start_node = 0;
    // variable that store the current node (at the beginning the current node is the start node)
    int current_node = start_node;
    
    int nind = 0;
    // indices of nodes except the current node
    for (int i = 0; i < n - cnt + 1; i++) if (i != current_node) indices[nind++] = i;
    
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
        (*objval) += distances[minim];
        
        // update the solution
        succ[current_node] = indices[minim];
        
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
        (*objval) += distances[minim];
        
        // update the solution
        succ[current_node] = indices[minim];
        
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
    
    // close the cycle with the last edge
    succ[current_node] = start_node;
    
    (*objval) += dist(current_node, start_node, inst);
    
    free(ind);
    free(distances);
    free(indices);
}
