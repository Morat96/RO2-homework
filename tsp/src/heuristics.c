//
//  heuristics.c
//  cplex
//
//  Created by Matteo Moratello on 04/05/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//

#include "heuristics.h"

// ************************** MODELS DEFINITION ************************* //
//
// - Nearest Neighbourhood (greedy)
// - Grasp (randomization)
// - Insertion
// - Insertion with Convex Hull
//
// ********************************************************************** //

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

// compute maximum of an array and return the corresponding index
int max(double *array, int arr_size) {
    
    int max_index = 0;
    double global_max = 0;
    
    if (arr_size == 1) return max_index;
    
    for (int i = 0; i < arr_size; i++) {
        if ( array[i] > global_max) {
            max_index = i;
            global_max = array[i];
        }
    }
    return max_index;
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
    free(global_min);
    return 0;
}

// reverse the direction of a part of the tour
void reverse_segment(instance* inst, int start, int end, int* succ) {
    
    int* index = (int*) calloc(inst -> nnodes, sizeof(int));
    int ind = succ[start];
    int cnt = 0;
    while (ind != end) {
        index[cnt++] = ind;
        ind = succ[ind];
    }
    if (cnt == 0) {
        succ[end] = start;
    }
    else {
        succ[index[0]] = start;
        for (int i = 1; i < cnt ; i++) {
            succ[index[i]] = index[i - 1];
        }
        succ[end] = index[cnt - 1];
    }
    free(index);
}

// order indices w.r.t. tour direction
void reorder(instance* inst, int f, int* s, int* t, int* succ) {
    
    int ind = f;
    for (int i = 0; i < inst -> nnodes; i++) {
        if (ind == *s) return;
        if (ind == *t) {
            *t = *s;
            *s = ind;
            return;
        }
        ind = succ[ind];
    }
}

// Refining algorithm
// 3-OPT move
void threeOpt(instance* inst, double* xstar) {
    
    printf("Refining Algorithm: 3-OPT move\n");
    double t1 = second();
    
    int *succ = (int *) calloc(inst -> nnodes, sizeof(int));
    int *comp = (int *) calloc(inst -> nnodes, sizeof(int));
    int ncomp = 0;
    
    build_sol(xstar, inst, succ, comp, &ncomp);
    
    while(1) {
        
        // first edge
        int nFirst = 0;
        int succFirst = 0;
        // second edge
        int nSecond = 0;
        int succSecond = 0;
        // third edge
        int nThird = 0;
        int succThird = 0;
        // case
        int move_case = 0;
        // objective function
        int delta = INT_MAX;
        int min_delta = INT_MAX;
        
        // Divide the tour in three segments deleting three different edges
        // and obtain a new tour combining the three segments
        // 7 ways to merge segments:
        // Case 1-2-3: merge a combination of only 2 segments (equivalent to two subsequent 2-opt moves)
        // Case 4-5-6: merge a combination of 3 segments (equivalent to three subsequent 2-opt moves)
        // Case 7: reverse the original tour merging all segments
        for (int i = 0; i < inst -> nnodes; i++) {
            for (int j = i + 1; j < inst -> nnodes; j++) {
                for (int z = j + 1; z < inst -> nnodes; z++) {
                    if (i != j && i != z && j != z) {
                        // reorder pairs of indices w.r.t. direction of the graph
                        reorder(inst, i, &j, &z, succ);
                        // A -> first tour segment , B -> second tour segment, C -> third tour segment
                        // First case: A'BC
                        delta = dist(z, i, inst) + dist(succ[z], succ[i], inst) - dist(i, succ[i], inst) - dist(z, succ[z], inst);
                        if (delta < min_delta) {
                            nFirst = i;
                            succFirst = succ[i];
                            nSecond = j;
                            succSecond = succ[j];
                            nThird = z;
                            succThird = succ[z];
                            move_case = 1;
                            min_delta = delta;
                        }
                        // Second case: ABC'
                        delta = dist(j, z, inst) + dist(succ[j], succ[z], inst) - dist(j, succ[j], inst) - dist(z, succ[z], inst);
                        if (delta < min_delta) {
                            nFirst = i;
                            succFirst = succ[i];
                            nSecond = j;
                            succSecond = succ[j];
                            nThird = z;
                            succThird = succ[z];
                            move_case = 2;
                            min_delta = delta;
                        }
                        // Third case: AB'C
                        delta = dist(i, j, inst) + dist(succ[i], succ[j], inst) - dist(i, succ[i], inst) - dist(j, succ[j], inst);
                        if (delta < min_delta) {
                            nFirst = i;
                            succFirst = succ[i];
                            nSecond = j;
                            succSecond = succ[j];
                            nThird = z;
                            succThird = succ[z];
                            move_case = 3;
                            min_delta = delta;
                        }
                        // Forth case: AB'C'
                        delta = dist(i, j, inst) + dist(succ[i], z, inst) + dist(succ[j], succ[z], inst) - dist(i, succ[i], inst) - dist(j, succ[j], inst) - dist(z, succ[z], inst);
                        if (delta < min_delta) {
                            nFirst = i;
                            succFirst = succ[i];
                            nSecond = j;
                            succSecond = succ[j];
                            nThird = z;
                            succThird = succ[z];
                            move_case = 4;
                            min_delta = delta;
                        }
                        // Fifth case: A'B'C
                        delta = dist(i, z, inst) + dist(succ[z], j, inst) + dist(succ[i], succ[j], inst) - dist(i, succ[i], inst) - dist(j, succ[j], inst) - dist(z, succ[z], inst);
                        if (delta < min_delta) {
                            nFirst = i;
                            succFirst = succ[i];
                            nSecond = j;
                            succSecond = succ[j];
                            nThird = z;
                            succThird = succ[z];
                            move_case = 5;
                            min_delta = delta;
                        }
                        // Sixth case: A'BC'
                        delta = dist(j, z, inst) + dist(succ[j], i, inst) + dist(succ[i], succ[z], inst) - dist(i, succ[i], inst) - dist(j, succ[j], inst) - dist(z, succ[z], inst);
                        if (delta < min_delta) {
                            nFirst = i;
                            succFirst = succ[i];
                            nSecond = j;
                            succSecond = succ[j];
                            nThird = z;
                            succThird = succ[z];
                            move_case = 6;
                            min_delta = delta;
                        }
                        // Seventh case: A'B'C'
                        delta = dist(i, succ[j], inst) + dist(succ[i], z, inst) + dist(j, succ[z], inst) - dist(i, succ[i], inst) - dist(j, succ[j], inst) - dist(z, succ[z], inst);
                        if (delta < min_delta) {
                            nFirst = i;
                            succFirst = succ[i];
                            nSecond = j;
                            succSecond = succ[j];
                            nThird = z;
                            succThird = succ[z];
                            move_case = 7;
                            min_delta = delta;
                        }
                    }
                }
            }
        }
        
        if (min_delta >= 0) break;
        
        // 3-OPT move based on best obj value
        switch (move_case) {
            case 1:
                // First case: A'BC
                reverse_segment(inst, succThird, nFirst, succ);
                succ[nThird] = nFirst;
                succ[succThird] = succFirst;
                break;
                
            case 2:
                // Second case: ABC'
                reverse_segment(inst, succSecond, nThird, succ);
                succ[nSecond] = nThird;
                succ[succSecond] = succThird;
                break;
                
            case 3:
                // Third case: AB'C
                reverse_segment(inst, succFirst, nSecond, succ);
                succ[nFirst] = nSecond;
                succ[succFirst] = succSecond;
                break;
            case 4:
                // Forth case: AB'C'
                reverse_segment(inst, succFirst, nSecond, succ);
                reverse_segment(inst, succSecond, nThird, succ);
                succ[nFirst] = nSecond;
                succ[succFirst] = nThird;
                succ[succSecond] = succThird;
                break;
            case 5:
                // Fifth case: A'B'C
                reverse_segment(inst, succThird, nFirst, succ);
                reverse_segment(inst, succFirst, nSecond, succ);
                succ[nThird] = nFirst;
                succ[succThird] = nSecond;
                succ[succFirst] = succSecond;
                break;
            case 6:
                // Sixth case: A'BC'
                reverse_segment(inst, succThird, nFirst, succ);
                reverse_segment(inst, succSecond, nThird, succ);
                succ[nSecond] = nThird;
                succ[succSecond] = nFirst;
                succ[succThird] = succFirst;
                break;
            case 7:
                // Seventh case: A'B'C'
                reverse_segment(inst, succThird, nFirst, succ);
                reverse_segment(inst, succFirst, nSecond, succ);
                reverse_segment(inst, succSecond, nThird, succ);
                succ[succSecond] = nFirst;
                succ[succThird] = nSecond;
                succ[succFirst] = nThird;
                break;
            default:
                break;
        }
    }
    
    double t2 = second();
    
    double objval = 0.0;
    for (int i = 0; i < inst -> nnodes; i ++) objval += dist(i, succ[i], inst);
    printf("Objective function value: %lf\n", objval);
    printf("3-OPT move time: %lf\n\n", t2 - t1);
    //print_solution_light(inst, succ);
    
    free(comp);
    free(succ);
}

// Refining algorithm
// 2-OPT move
void twOpt(instance* inst, double* xstar) {
    
    printf("Refining Algorithm: 2-OPT move\n");
    double t1 = second();
    
    int *succ = (int *) calloc(inst -> nnodes, sizeof(int));
    int *comp = (int *) calloc(inst -> nnodes, sizeof(int));
    int ncomp = 0;
    
    build_sol(xstar, inst, succ, comp, &ncomp);
    
    int *index = (int *) calloc(inst -> nnodes, sizeof(int));
    
    while(1) {
        
        // first edge
        int nFirst = 0;
        int succ_f = 0;
        // second edge
        int nSecond = 0;
        int succ_s = 0;
        // objective function
        int delta = INT_MAX;
        int min_delta = INT_MAX;
        
        // for each couple of edges, compute obj function value and pick the minimim one
        // store indices of edges with the lowest obj func value
        for (int i = 0; i < inst -> nnodes; i++) {
            for (int j = i + 1; j < inst -> nnodes; j++) {
                delta = dist(i, j, inst) + dist(succ[i], succ[j], inst) - dist(i, succ[i], inst) - dist(j, succ[j], inst);
                if (delta < min_delta) {
                    nFirst = i;
                    succ_f = succ[i];
                    nSecond = j;
                    succ_s = succ[j];
                    min_delta = delta;
                }
            }
        }
        
        // if no improvement, exit
        if (min_delta >= 0) break;
        
        // change link of edges and their orientations w.r.t. optimal obj. func. value
        int ind = succ[succ_f];
        int cnt = 0;
        while (ind != nSecond) {
            index[cnt++] = ind;
            ind = succ[ind];
        }
        if (cnt == 0) {
            succ[ind] = succ_s;
            succ[nSecond] = succ_f;
        }
        else {
            succ[index[0]] = succ_f;
            for (int i = 1; i < cnt ; i++) {
                succ[index[i]] = index[i - 1];
            }
            succ[nSecond] = index[cnt - 1];
        }
        // crossing of straight lines
        // (a,a'), (b,b') -> (a,b), (a',b')
        // where (a,a') is the first edge and (b,b') is the second one.
        succ[nFirst] = nSecond;
        succ[succ_f] = succ_s;
    }
    
    double t2 = second();
    
    double objval = 0.0;
    for (int i = 0; i < inst -> nnodes; i ++) objval += dist(i, succ[i], inst);
    printf("Objective function value: %lf\n", objval);
    printf("2-OPT move time: %lf\n\n", t2 - t1);
    //print_solution_light(inst, succ);
    
    free(index);
    free(succ);
    free(comp);
}

// Insertion heuristic with Convex Hull
void insertion_ch(instance *inst, double *xstar) {
    
    printf("Resolve instance \"%s\" with Insertion with Convex Hull\n\n", inst -> input_file);
    double t1 = second();
    
    // number of nodes
    int n = inst -> nnodes;
    
    int *indices = (int *) calloc(n - 1, sizeof(int));
    
    // solution nodes
    int *sol = (int *) calloc(n, sizeof(int));
    
    // number of nodes in solution
    int n_nodes_sol = 0;
    
    // save all nodes in an array of Point
    Point *points = (Point *) calloc(n, sizeof(Point));
    for (int i = 0; i < n; i++) {
        points[i].x = inst -> xcoord[i];
        points[i].y = inst -> ycoord[i];
    }
    
    Point *ch = (Point *) calloc(n, sizeof(Point));
    
    // size of the convex hull
    int size = 0;
    
    // compute the convex hull of the points
    // O(nLogn) complexity
    if (convexHull(points, n, ch, &size)) print_error("error in build the Convex Hull");
    
    // save the current solution
    int nnodes = 0;
    for (int j = 0; j < size; j++) {
        for (int i = 0; i < n; i++) {
            if (inst -> xcoord[i] == ch[j].x && inst -> ycoord[i] == ch[j].y) sol[nnodes++] = i;
        }
    }
    
    // indices of all other nodes
    int counter = 0;
    int n_ind = 0;
    for (int i = 0; i < n; i++) {
        counter = 0;
        for (int j = 0; j < size; j++) {
            if (inst -> xcoord[i] == ch[j].x && inst -> ycoord[i] == ch[j].y) counter ++;
        }
        if (!counter) indices[n_ind++] = i;
    }
    
    double t_ch = second();
    printf("Convex Hull build time: %lf\n", t_ch - t1);
    
    int start = 0;
    n_nodes_sol = size;
    
    // End first part
    // Second part
    // insert nodes on the solution until all edges belong to C
    
    int cnt = size;
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
                double cost_h = dist(sol[c], indices[i], inst) + dist(indices[i], sol[c + 1], inst) - dist(sol[c], sol[c + 1], inst);
                
                // updating the best cost value
                if (cost_h < best_cost) {
                    best_h = indices[i];
                    best_h_index = i;
                    best_pos = c;
                    best_cost = cost_h;
                }
            }
        }
        
        // last edge (the one that closes the cycle)
        for (int i = 0; i < n - cnt; i++) {
            
            double cost_h = dist(sol[n_nodes_sol - 1], indices[i], inst) + dist(indices[i], sol[start], inst) - dist(sol[n_nodes_sol - 1], sol[start], inst);
            
            // updating the best cost value
            if (cost_h < best_cost) {
                best_h = indices[i];
                best_h_index = i;
                best_pos = n_nodes_sol - 1;
                best_cost = cost_h;
            }
        }
        
        // update the solution with the best edge
        for (int c = n - 1; c > best_pos; c--) sol[c] = sol[c - 1];
        sol[best_pos + 1] = best_h;
        
        // remove the previous node
        for (int c = best_h_index; c < n - cnt; c++) indices[c] = indices[c + 1];
        
        cnt ++;
        n_nodes_sol ++;
    }
    
    // build the xstar and the objective function value
    for (int i = 0; i < n - 1; i++) {
        objval += dist(sol[i], sol[i+1], inst);
        xstar[xpos(sol[i], sol[i+1], inst)] = 1.0;
    }
    
    xstar[xpos(sol[n - 1], sol[start], inst)] = 1.0;
    objval += dist(sol[n - 1], sol[start], inst);
    
    printf("Objective function value: %lf\n", objval);
    //for ( int j = 0; j < ncols; j++ ) printf(" ... qstar[%3d] = %10.2lf \n", j+1, xstar[j]);
    
    // build and print the solution
    int *succ = (int *) calloc(inst->nnodes, sizeof(int));
    int *comp = (int *) calloc(inst->nnodes, sizeof(int));
    int ncomp = 0;
    
    build_sol(xstar, inst, succ, comp, &ncomp);
    
    double t2 = second();
    printf("Insertion CH time: %lf\n\n", t2 - t1);
    
    //print_solution_light(inst, succ);
    
    free(sol);
    free(indices);
    free(ch);
    free(points);
    
    free(comp);
    free(succ);
}

// Insertion heuristic
void insertion(instance *inst, double *xstar) {
    
    printf("Resolve instance \"%s\" with Insertion\n\n", inst -> input_file);
    double t1 = second();
    
    // number of nodes
    int n = inst -> nnodes;
    
    int *indices = (int *) calloc(n - 1, sizeof(int));
    double *distances = (double *) calloc(n - 1, sizeof(double));
    
    // solution nodes
    int *sol = (int *) calloc(n, sizeof(int));
    
    // number of nodes in solution
    int n_nodes_sol = 0;
    
    // start node
    int start_node = 0;
    
    // index of first node
    int start = 0;
    
    sol[n_nodes_sol++] = start_node;
    
    // compute distances between start node and all the other nodes
    int nind = 0;
    for (int i = 0; i < n; i++) if (i != start_node) indices[nind++] = i;
    for (int i = 0; i < n - 1; i++) distances[i] = dist(start, indices[i], inst);
    
    // search the farthest node w.r.t. start node
    int maxim = max(distances, n - 1);
    
    // second node in the solution
    sol[n_nodes_sol++] = indices[maxim];
    
    // remove the previous node
    for (int c = maxim; c < n - 1; c++) indices[c] = indices[c+1];
    
    // End first part
    // Second part
    // insert nodes on the solution until all edges belong to C
    
    int cnt = 2;
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
                double cost_h = dist(sol[c], indices[i], inst) + dist(indices[i], sol[c + 1], inst) - dist(sol[c], sol[c + 1], inst);
                
                // updating the best cost value
                if (cost_h < best_cost) {
                    best_h = indices[i];
                    best_h_index = i;
                    best_pos = c;
                    best_cost = cost_h;
                }
            }
        }
        
        // last edge (the one that closes the cycle)
        for (int i = 0; i < n - cnt; i++) {
            
            double cost_h = dist(sol[n_nodes_sol - 1], indices[i], inst) + dist(indices[i], sol[start], inst) - dist(sol[n_nodes_sol - 1], sol[start], inst);
            
            // updating the best cost value
            if (cost_h < best_cost) {
                best_h = indices[i];
                best_h_index = i;
                best_pos = n_nodes_sol - 1;
                best_cost = cost_h;
            }
        }
        
        // update the solution with the best edge
        for (int c = n - 1; c > best_pos; c--) sol[c] = sol[c - 1];
        sol[best_pos + 1] = best_h;
        
        // remove the previous node
        for (int c = best_h_index; c < n - cnt; c++) indices[c] = indices[c + 1];
        
        cnt ++;
        n_nodes_sol ++;
    }
    
    // build the xstar and the objective function value
    for (int i = 0; i < n - 1; i++) {
        objval += dist(sol[i], sol[i+1], inst);
        xstar[xpos(sol[i], sol[i+1], inst)] = 1.0;
    }
    
    xstar[xpos(sol[n - 1], sol[start], inst)] = 1.0;
    objval += dist(sol[n - 1], sol[start], inst);
    
    printf("Objective function value: %lf\n", objval);
    //for ( int j = 0; j < ncols; j++ ) printf(" ... qstar[%3d] = %10.2lf \n", j+1, xstar[j]);
    
    // build and print the solution
    int *succ = (int *) calloc(inst->nnodes, sizeof(int));
    int *comp = (int *) calloc(inst->nnodes, sizeof(int));
    int ncomp = 0;
    
    build_sol(xstar, inst, succ, comp, &ncomp);
    
    double t2 = second();
    printf("Insertion time: %lf\n\n", t2 - t1);
    
    //print_solution_light(inst, succ);
    
    free(sol);
    free(indices);
    free(distances);
    
    free(comp);
    free(succ);
}

// Nearest Neighborhood (greedy)
void NearNeigh(instance *inst, double *xstar) {
    
    printf("Resolve instance \"%s\" with Nearest Neighborhood\n\n", inst -> input_file);
    double t1 = second();
    
    // number of variables
    
    // number of nodes of the instance
    int n = inst -> nnodes;
    
    // solution
    //double *xstar = (double *) calloc(ncols, sizeof(double));
    
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
    
    printf("Objective function value: %lf\n", objval);
    //for ( int j = 0; j < ncols; j++ ) printf(" ... qstar[%3d] = %10.2lf \n", j+1, xstar[j]);
    
    free(distances);
    free(indices);
    
    // build and print the solution
    int *succ = (int *) calloc(inst->nnodes, sizeof(int));
    int *comp = (int *) calloc(inst->nnodes, sizeof(int));
    int ncomp = 0;
    
    build_sol(xstar, inst, succ, comp, &ncomp);
    
    double t2 = second();
    printf("Nearest Neighborhood time: %lf\n\n", t2 - t1);
    
    //print_solution(inst, succ);
    
    free(comp);
    free(succ);
}

// Grasp (randomization)
void grasp(instance *inst, double *xstar) {
    
    printf("Resolve instance \"%s\" with Grasp\n\n", inst -> input_file);
    double t1 = second();
    
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
    
    printf("Objective function value: %lf\n", objval);
    //for ( int j = 0; j < ncols; j++ ) printf(" ... qstar[%3d] = %10.2lf \n", j+1, xstar[j]);
    
    free(ind);
    free(distances);
    free(indices);
    
    // build and print the solution
    int *succ = (int *) calloc(inst->nnodes, sizeof(int));
    int *comp = (int *) calloc(inst->nnodes, sizeof(int));
    int ncomp = 0;
    
    build_sol(xstar, inst, succ, comp, &ncomp);
    
    double t2 = second();
    printf("Grasp time: %lf\n\n", t2 - t1);
    
    //print_solution(inst, succ);
    
    free(comp);
    free(succ);
    
}

// Grasp cycled (not a good idea for me)
void grasp_cycled(instance *inst) {
    
    printf("Resolve instance \"%s\" with Grasp\n\n", inst -> input_file);
    
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
