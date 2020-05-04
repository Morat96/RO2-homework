//
//  heuristics.c
//  cplex
//
//  Created by Matteo Moratello on 04/05/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//

#include "heuristics.h"

double dist(int i, int j, instance *inst);
int xpos(int i, int j, instance *inst);
void build_sol(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp);
void print_solution(instance *inst, int *succ);
void print_solution_light(instance *inst, int *succ);
void print_error(const char *err);

// compute minimum of an array and return the corresponding index
int min(double *array, int arr_size) {
    
    int min_index = 0;
    double global_min = INT_MAX;
    
    if (arr_size == 1) return min_index;
    
    for (int i = 0; i < arr_size; i++) {
        if ( array[i] < global_min) {
            min_index = i;
            global_min = array[i];
        }
    }
    return min_index;
}

// compute the three smaller values of an array and store the indexes in the "ind" array
// return 1 if the array size is smaller than three, otherwise return 0
int three_min(double *array, int arr_size, int *ind) {
    
    // There should be at least three elements
    if (arr_size < 3) {
        printf(" Invalid Input ");
        return 1;
    }
    
    double *global_min = (double *) calloc(3, sizeof(double));
    for (int i = 0; i < 3; i++) global_min[i] = INT_MAX;
    
    // array structure
    // ind[0] = global minimum
    // ind[1] = second minimum
    // ind[2] = third minimum
    for (int i = 0; i < arr_size; i++) {
        if (array[i] < global_min[0]) {
            global_min[2] = global_min[1];
            global_min[1] = global_min[0];
            global_min[0] = array[i];
            ind[2] = ind[1];
            ind[1] = ind[0];
            ind[0] = i;
        }
        else if (array[i] < global_min[1]) {
            global_min[2] = global_min[1];
            global_min[1] = array[i];
            ind[2] = ind[1];
            ind[1] = i;
        }
        else if (array[i] < global_min[2]) {
            global_min[2] = array[i];
            ind[2] = i;
        }
    }
    return 0;
}

// Nearest Neighborhood (greedy)
void NearNeigh(instance *inst) {
    
    printf("Resolve instance \"%s\" with Nearest Neighborhood\n\n", inst -> input_file);
    
    // number of variables
    int ncols = ((inst -> nnodes)*(inst -> nnodes - 1)) / 2;
    // number of nodes of the instance
    int n = inst -> nnodes;
    
    // solution
    double *xstar = (double *) calloc(ncols, sizeof(double));
    
    // counter of resolved nodes
    int cnt = 1;
    
    int *indices = (int *) calloc(inst -> nnodes - 1, sizeof(int));
    double *distances = (double *) calloc(inst -> nnodes - 1, sizeof(double));
    
    // first node [0, nnodes - 1]
    int start_node = 0;
    // variable that store the current node (at the beginning the current node is the start node)
    int current_node = start_node;
    
    int nind = 0;
    // indices of nodes except the current node
    for (int i = 0; i < n - cnt + 1; i++) if (i != current_node) indices[nind++] = i;
    
    // obj function value
    double objval = 0.0;
    
    // in each cycle find the minimum edge
    while (cnt < inst -> nnodes) {
    
        // compute distances
        for (int i = 0; i < n - cnt; i++) distances[i] = dist(current_node, indices[i], inst);
        
        // find the shorter edge from the current node
        int minim = min(distances, n - cnt);
        
        // update objective function value
        objval += distances[minim];
        
        // update the solution
        xstar[xpos(current_node, indices[minim], inst)] = 1.0;
        
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
    xstar[xpos(current_node, start_node, inst)] = 1.0;
    objval += dist(current_node, start_node, inst);
    
    printf("Objective function value: %lf\n" , objval);
    //for ( int j = 0; j < ncols; j++ ) printf(" ... qstar[%3d] = %10.2lf \n", j+1, xstar[j]);
    
    free(distances);
    free(indices);

    // build and print the solution
    int *succ = (int *) calloc(inst->nnodes, sizeof(int));
    int *comp = (int *) calloc(inst->nnodes, sizeof(int));
    int ncomp = 0;
    
    build_sol(xstar, inst, succ, comp, &ncomp);
    print_solution(inst, succ);
    
    free(xstar);
    free(comp);
    free(succ);
}

// Grasp (randomization)
void grasp(instance *inst) {
    
    printf("Resolve instance \"%s\" with GRASP\n\n", inst -> input_file);
    
    // number of variables
    int ncols = ((inst -> nnodes)*(inst -> nnodes - 1)) / 2;
    // number of nodes of the instance
    int n = inst -> nnodes;
    
    // solution
    double *xstar = (double *) calloc(ncols, sizeof(double));
    
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
    
    // obj function value
    double objval = 0.0;
    
    // randomize seed
    srand((unsigned int)time(NULL));
    
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
        if (rv <= 0.8) minim = ind[0];
        else if (rv <= 0.9) minim = ind[1];
        else minim = ind[2];
        
        // update objective function value
        objval += distances[minim];
        
        // update the solution
        xstar[xpos(current_node, indices[minim], inst)] = 1.0;
        
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
        xstar[xpos(current_node, indices[minim], inst)] = 1.0;
        
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
    xstar[xpos(current_node, start_node, inst)] = 1.0;
    objval += dist(current_node, start_node, inst);
    
    printf("Objective function value: %lf\n" , objval);
    //for ( int j = 0; j < ncols; j++ ) printf(" ... qstar[%3d] = %10.2lf \n", j+1, xstar[j]);
    
    free(ind);
    free(distances);
    free(indices);
    
    // build and print the solution
    int *succ = (int *) calloc(inst->nnodes, sizeof(int));
    int *comp = (int *) calloc(inst->nnodes, sizeof(int));
    int ncomp = 0;
    
    build_sol(xstar, inst, succ, comp, &ncomp);
    print_solution(inst, succ);
    
    free(xstar);
    free(comp);
    free(succ);
    
}

// Grasp cycled (not a good idea for me)
void grasp_cycled(instance *inst) {
    
    printf("Resolve instance \"%s\" with GRASP\n\n", inst -> input_file);
    
    // number of variables
    int ncols = ((inst -> nnodes)*(inst -> nnodes - 1)) / 2;
    // number of nodes of the instance
    int n = inst -> nnodes;
    
    // solution
    double *xstar = (double *) calloc(ncols, sizeof(double));
    double *best_x_star = (double *) calloc(ncols, sizeof(double));
    double best_objval = CPX_INFBOUND;
    
    int *indices = (int *) calloc(inst -> nnodes - 1, sizeof(int));
    double *distances = (double *) calloc(inst -> nnodes - 1, sizeof(double));
    int *ind = (int *) calloc(3, sizeof(int));
    
    // randomize seed
    srand((unsigned int)time(NULL));
    
    for (int cyc = 0; cyc < n; cyc ++) {
        
        // counter of resolved nodes
        int cnt = 1;
        
        // first node [0, nnodes - 1]
        int start_node = cyc;
        // variable that store the current node (at the beginning the current node is the start node)
        int current_node = start_node;
        
        int nind = 0;
        // indices of nodes except the current node
        for (int i = 0; i < n; i++) if (i != current_node) indices[nind++] = i;
        
        // obj function value
        for (int x = 0; x < ncols; x++) xstar[x] = 0.0;
        double objval = 0.0;
        
        // in each cycle find the minimum edge
        while (cnt < inst -> nnodes - 2) {
            
            // compute distances
            for (int i = 0; i < n - cnt; i++) distances[i] = dist(current_node, indices[i], inst);
            
            // find the shorter edge from the current node
            int minim = 0;
            if (three_min(distances, n - cnt, ind)) print_error("Error in function three_min");
            
            /* choose the next node with weighted probability.
             The node corresponding to the smaller edge has 50% probability to be chosen.
             The other two nodes have 25% probability to be chosen. */
            
            float rv = ((float)rand()/(float)(RAND_MAX));
            
            if (rv <= 0.6) minim = ind[0];
            else if (rv <= 0.8) minim = ind[1];
            else minim = ind[2];
            
            // update objective function value
            objval += distances[minim];
            
            // update the solution
            xstar[xpos(current_node, indices[minim], inst)] = 1.0;
            
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
        for (int fin = 0; fin < 2; fin++) {
            
            // compute distances
            for (int i = 0; i < n - cnt; i++) distances[i] = dist(current_node, indices[i], inst);
            
            // find the shorter edge from the current node
            int minim = min(distances, n - cnt);
            
            // update objective function value
            objval += distances[minim];
            
            // update the solution
            xstar[xpos(current_node, indices[minim], inst)] = 1.0;
            
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
        xstar[xpos(current_node, start_node, inst)] = 1.0;
        objval += dist(current_node, start_node, inst);
        
        // update objective function value and xstar
        if (objval < best_objval) {
            best_objval = objval;
            for (int x = 0; x < ncols; x++) best_x_star[x] = xstar[x];
        }
    }
    
    free(ind);
    free(distances);
    free(indices);
    
    // build and print the solution
    int *succ = (int *) calloc(inst->nnodes, sizeof(int));
    int *comp = (int *) calloc(inst->nnodes, sizeof(int));
    int ncomp = 0;
    
    printf("Best Objective function value: %lf\n" , best_objval);
    build_sol(best_x_star, inst, succ, comp, &ncomp);
    print_solution(inst, succ);
    
    free(xstar);
    free(comp);
    free(succ);
    
}
