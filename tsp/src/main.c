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
void twOpt(instance* inst, double* xstar);

void print_error(const char *err) { printf("\n\n ERROR: %s \n\n", err); fflush(NULL); exit(1); }


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
    twOpt(&inst, xstar);
    
    //if ( TSPopt(&inst, t1) ) print_error(" error within TSPopt()");
    double t2 = second();
    
    if ( VERBOSE >= 1 ) printf("\n... TSP problem solved in %lf sec\n", t2 - t1);
    
    free_instance(&inst);
    
    return 0;
}



