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
void print_solution_light(instance *inst, int *succ);
void print_error(const char *err) { printf("\n\n ERROR: %s \n\n", err); fflush(NULL); exit(1); }
int xpos(int i, int j, instance *inst);
double dist(int i, int j, instance *inst);

// others
void smallerKnodes(instance* inst, int** distances);
void sort(double* array, int* ind, int begin, int end);
void swap(double* first, double* second);
void swap_int(int* first, int* second);

// metaheuristics
void vns(instance* inst, int iter, int k);
void tabu_search(instance* inst, int iter, int list_size);
void simulated_annealing(instance* inst, int iter, int size);
void genetic_algorithm(instance* inst, int generations, int population_size, int crossover_size);
void multi_start(instance* inst, int iter);

/*void NearNeigh(instance *inst, double *xstar);
void grasp(instance *inst, double *xstar);
void insertion(instance *inst, double *xstar);
void insertion_ch(instance *inst, double *xstar);
void random_solution(instance* inst, double* xstar);
void twOpt(instance* inst, double* xstar);
void threeOpt(instance* inst, double* xstar);*/

/**
 Free the struct instance for TSP.

 @param inst instance of the struct "instance" for TSP problem.
 */
void free_instance(instance *inst)
{
    free(inst -> xcoord);
    free(inst -> ycoord);
    for (int i = 0; i < 4; i++) free(inst -> sol_thread[i]);
    free(inst -> sol_thread);
}

/**
 Main.

 @param argc number of arguments.
 @param argv array of char with arguments.
 @return 0 if succeeded, 1 otherwise.
 */
int main(int argc, char **argv) {
    
    if ( argc < 2 ) { printf("Usage: %s -help for help\n", argv[0]); exit(1); }
    if ( VERBOSE >= 2 ) { for (int a = 0; a < argc; a++) printf("%s ", argv[a]); printf("\n"); }
    
    double t1 = second();
    instance inst;
    
    parse_command_line(argc, argv, &inst);
    read_input(&inst);
    
    int ncols = ((inst.nnodes) * (inst.nnodes - 1)) / 2;
    double *xstar = (double *) calloc(ncols, sizeof(double));
    
    //NearNeigh(&inst, xstar);
    //grasp(&inst, xstar);
    //insertion(&inst, xstar);
    //insertion_ch(&inst, xstar);
    //random_solution(&inst, xstar);
    //twOpt(&inst, xstar);
    //threeOpt(&inst, xstar);
    //smallerKnodes(&inst);
    
    multi_start(&inst, 10);
    if (inst.vns) vns(&inst, 100, 5);
    else if (inst.tabu_search) tabu_search(&inst, 5000, 100);
    else if (inst.sim_annealing) simulated_annealing(&inst, 500, 50);
    else if (inst.genetic) genetic_algorithm(&inst, 2, 100, 20);
    else {
        //if ( TSPopt(&inst, t1) ) print_error(" error within TSPopt()");
    }
    
    double t2 = second();
    
    if ( VERBOSE >= 1 ) printf("\n... TSP problem solved in %lf sec\n", t2 - t1);
    
    free(xstar);
    free_instance(&inst);
    
    return 0;
}

//////////////////////
/////// Others ///////
//////////////////////
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
