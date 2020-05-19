//
//  main.c
//  cplex
//
//  Created by Matteo Moratello on 12/03/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//
#include "tsp.h"

double second(void);
void read_input(instance *inst);
void parse_command_line(int argc, char** argv, instance *inst);
int TSPopt(instance *inst, double t1);
void NearNeigh(instance *inst, double *xstar);
void grasp(instance *inst, double *xstar);
void insertion(instance *inst, double *xstar);
void insertion_ch(instance *inst, double *xstar);
void random_solution(instance* inst, double* xstar);
void twOpt(instance* inst, double* xstar);
void threeOpt(instance* inst, double* xstar);
void print_solution_light(instance *inst, int *succ);
void print_error(const char *err) { printf("\n\n ERROR: %s \n\n", err); fflush(NULL); exit(1); }
int xpos(int i, int j, instance *inst);
double dist(int i, int j, instance *inst);
void smallerKnodes(instance* inst, int** distances);
void sort(double* array, int* ind, int begin, int end);
void swap(double* first, double* second);
void swap_int(int* first, int* second);

void free_instance(instance *inst)
{
    free(inst->xcoord);
    free(inst->ycoord);
    for (int i = 0; i < 4; i++) free(inst -> sol_thread[i]);
    free(inst->sol_thread);
}

int main(int argc, char **argv) {
    
    if ( argc < 2 ) { printf("Usage: %s -help for help\n", argv[0]); exit(1); }
    if ( VERBOSE >= 2 ) { for (int a = 0; a < argc; a++) printf("%s ", argv[a]); printf("\n"); }
    
    double t1 = second();
    instance inst;
    
    parse_command_line(argc, argv, &inst);
    read_input(&inst);
    
    int ncols = ((inst.nnodes) * (inst.nnodes - 1)) / 2;
    double *xstar = (double *) calloc(ncols, sizeof(double));
    
    NearNeigh(&inst, xstar);
    //grasp(&inst, xstar);
    //insertion(&inst, xstar);
    //insertion_ch(&inst, xstar);
    //random_solution(&inst, xstar);
    twOpt(&inst, xstar);
    threeOpt(&inst, xstar);
    //smallerKnodes(&inst);
    
    //if ( TSPopt(&inst, t1) ) print_error(" error within TSPopt()");
    double t2 = second();
    
    if ( VERBOSE >= 1 ) printf("\n... TSP problem solved in %lf sec\n", t2 - t1);
    
    free(xstar);
    free_instance(&inst);
    
    return 0;
}

// build a TSP random solution
void random_solution(instance* inst, double* xstar) {
    
    printf("Resolve instance \"%s\" with Random\n\n", inst -> input_file);
    
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

    double objval = 0.0;
    
    for (int i = 0; i < inst -> nnodes - 1; i++) {
        xstar[xpos(sol[i], sol[i + 1], inst)] = 1.0;
        objval += dist(sol[i], sol[i + 1], inst);
    }
    xstar[xpos(sol[inst -> nnodes - 1], sol[0], inst)] = 1.0;
    objval += dist(sol[inst -> nnodes - 1], sol[0], inst);
    
    printf("Objective function value: %lf\n\n", objval);
    
    free(sol);
    free(index);
}

void smallerKnodes(instance* inst, int** distances) {
    
    int cnt = 1;

    cnt = 1;
    double* dis = (double*) calloc(inst -> nnodes - 1, sizeof(double));
    
    for (int i = 0; i < inst -> nnodes - 1; i++) {
        int c = 0;
        for (int j = i + 1; j < inst -> nnodes; j++) {
            dis[i] = dist(i, j, inst);
            distances[i][c++] = j;
        }
        sort(dis, distances[i], 0, inst -> nnodes - cnt - 1);
        cnt ++;
    }
    free(dis);
}

// Quicksort
void sort(double* array, int* ind, int begin, int end) {
    int pivot, l, r;
    if (end > begin) {
        pivot = array[begin];
        l = begin + 1;
        r = end+1;
        while(l < r)
            if (array[l] < pivot)
                l++;
            else {
                r--;
                swap(&array[l], &array[r]);
                swap_int(&ind[l], &ind[r]);
            }
        l--;
        swap(&array[begin], &array[l]);
        swap_int(&ind[begin], &ind[l]);
        sort(array, ind, begin, l);
        sort(array, ind, r, end);
    }
}

void swap(double* first, double* second) {
    double temp = *first;
    *first = *second;
    *second = temp;
}

void swap_int(int* first, int* second) {
    int temp = *first;
    *first = *second;
    *second = temp;
}
