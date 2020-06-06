//
//  tabu_search.c
//  cplex
//
//  Created by Matteo Moratello on 04/06/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//

#include "tabu_search.h"

// *********************************** Tabù Search ********************************* //
//
// Description: TO DO.
//
// ********************************************************************************* //

/**
 Tabù search algorithm for TSP.

 @param inst instance of the struct "instance" for TSP problem.
 @param iter number of iterations.
 @param list_size tabù list size.
 */
void tabu_search(instance* inst, int iter, int list_size) {
    
    double best_obj = 0.0;
    double current_obj = 0.0;
    // incumbent
    int *succ_best = (int *) calloc(inst -> nnodes, sizeof(int));
    double* xstar = (double*) calloc(((inst -> nnodes * (inst -> nnodes - 1))/2), sizeof(double));
    
    // Phase 1
    // find a solution with NN and find a local optima using 2-OPT
    NearNeigh(inst, xstar);
    
    int *succ = (int *) calloc(inst -> nnodes, sizeof(int));
    int *comp = (int *) calloc(inst -> nnodes, sizeof(int));
    int ncomp = 0;
    
    build_sol(xstar, inst, succ, comp, &ncomp);
    
    // search in a neighborhood with 2-OPT algorithm
    twOpt(inst, succ, &current_obj);
    
    best_obj = current_obj;
    // update incumbent
    for (int i = 0; i < inst -> nnodes; i++) succ_best[i] = succ[i];
    
    // Phase 2
    // Tabù search
    tabu_vert* tabu_list = (tabu_vert*) calloc(list_size, sizeof(tabu_vert));
    
    int current_tabu_size = 0;
    int check = 0;
    int max_size_reached = 0;
    int n_empty_list = 0.8 * list_size;
    
    // algorithm iterations
    for (int i = 0; i < iter; i++) {
        findBestEdgesSwap(inst, succ, &current_obj, tabu_list, &current_tabu_size, list_size, check);
        if (current_obj < best_obj) {
            best_obj = current_obj;
            for (int i = 0; i < inst -> nnodes; i++) succ_best[i] = succ[i];
        }
        if (current_tabu_size == list_size) {
            current_tabu_size = 0;
            check = 1;
            max_size_reached = 1;
        }
        if (max_size_reached && (i % n_empty_list == 0)) {
            check = 0;
            current_tabu_size = 0;
            twOpt(inst, succ, &current_obj);
            if (current_obj < best_obj) {
                best_obj = current_obj;
                for (int i = 0; i < inst -> nnodes; i++) succ_best[i] = succ[i];
            }
        }
    }
    
    printf("\nBest objective function value: %lf\n", best_obj);
    free(tabu_list);
    free(succ);
    free(succ_best);
    free(comp);
    free(xstar);
}

/**
 Find the less worse edge swap in the TSP tour.

 @param inst instance of the struct "instance" for TSP problem.
 @param succ TPS solution with successors format.
 @param objval TSP objective function value.
 @param tabu_list tabù list.
 @param tabu_size current tabù list size.
 @param list_size original tabù list size.
 @param check if empty the list or not.
 */
void findBestEdgesSwap(instance* inst, int* succ, double* objval, tabu_vert* tabu_list, int* tabu_size, int list_size, int check) {
    
    int ls = 0;
    if (check) ls = list_size;
    else ls = (*tabu_size);
    (*objval) = 0.0;
    
    int *index = (int *) calloc(inst -> nnodes, sizeof(int));
    // first edge
    int nFirst = 0;
    int succ_f = 0;
    // second edge
    int nSecond = 0;
    int succ_s = 0;
    // objective function
    int delta = INT_MAX;
    int min_delta = INT_MAX;
    int jump = 0;
    
    // for each couple of edges, compute obj function value and pick the minimim one
    // store indices of edges with the lowest obj func value
    for (int i = 0, a = 0; a < inst -> nnodes - 1; a++, i = succ[i]) {
        for (int j = succ[i], b = a + 1; b < inst -> nnodes; b++, j = succ[j]) {
            for (int z = 0; z < ls; z++) {
                tabu_vert verts = tabu_list[z];
                if ((verts.vert1 == i || verts.vert2 == j)) {
                    jump = 1;
                    break;
                }
                jump = 0;
            }
            if (!jump) {
                delta = dist(i, j, inst) + dist(succ[i], succ[j], inst) - dist(i, succ[i], inst) - dist(j, succ[j], inst);
                if (delta != 0 && delta < min_delta) {
                    nFirst = i;
                    succ_f = succ[i];
                    nSecond = j;
                    succ_s = succ[j];
                    min_delta = delta;
                }
            }
        }
    }
    
    //printf("delta: %d\n", min_delta);
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
    
    tabu_vert verts;
    verts.vert1 = nFirst;
    verts.vert2 = nSecond;
    tabu_list[(*tabu_size)] = verts;
    (*tabu_size) ++;
    
    //double objval = 0.0;
    for (int i = 0; i < inst -> nnodes; i ++) (*objval) += dist(i, succ[i], inst);
    //printf("Objective function value: %lf\n", (*objval));
    
    free(index);
}
