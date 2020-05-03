//
//  loop_method.c
//  cplex
//
//  Created by Matteo Moratello on 13/04/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//

#include "loop_method.h"

double second(void);
void print_error(const char *err);
void print_solution(instance *inst, int *succ);
void build_sol(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp);
int xpos(int i, int j, instance *inst);

// ********************************************************************** //
// ***************************** LOOP METHOD **************************** //
// ********************************************************************** //


// ************************** MODEL DEFINITION ************************** //
//
// obj: compute an optimal solution adding SEC constraints at the end of each run of CPXmipopt
// First version
// Add SEC constraints until ncomp (number of components) is equal to 1
//
// Second version
// In a first phase set CPLEX parameters on order to return a non-optimal solution but
// use it in order to produce valid constraints for the problem. Once we obtain 1 component
// come back to the original CPLEX setting in order to obtain the optimal solution for the problem.
// Parameters to tune:
// CPXsetintparam(env, CPX_PARAM_INTSOLLIM, 3); --> incumbent
// CPXsetintparam(env, CPX_PARAM_EPGAP, 0.03);  --> solution gap
// CPXsetintparam(env, CPX_PARAM_NODELIM, 0);   --> limit on the branching node
//
// ********************************************************************** //
void loop_method(instance *inst, CPXENVptr env, CPXLPptr lp, double t1) {
    
    int *succ = (int *) calloc(inst->nnodes, sizeof(int));
    int *comp = (int *) calloc(inst->nnodes, sizeof(int));
    // number of components
    int ncomp = 9999;
    // iteration number
    int n = 1;
    // objective function value
    double objval = 0;
    double t2 = 0;
    // phases
    int first_phase = 1;
    int second_phase = 0;
    // set parameters
    //CPXsetintparam(env, CPX_PARAM_INTSOLLIM, 3);
    //CPXsetintparam(env, CPX_PARAM_EPGAP, 0.03);
    CPXsetintparam(env, CPX_PARAM_NODELIM, 0);
    
    printf("Resolve instance \"%s\" with Loop Method\n", inst -> input_file);
    printf("--------------------------------------\n");
    printf("             First Phase\n");
    printf("--------------------------------------\n");
    
    // Loop method (Benders)
    while (ncomp != 1) {
        
        // if phase 2 set default parameters
        if (second_phase == 1) {
            printf("--------------------------------------\n");
            printf("             Second Phase\n");
            printf("--------------------------------------\n");
            //CPXsetintparam(env, CPX_PARAM_INTSOLLIM, 2100000000);
            //CPXsetintparam(env, CPX_PARAM_EPGAP, 1);
            CPXsetintparam(env, CPX_PARAM_NODELIM, 2100000000);
            second_phase = 0;
        }
        
        printf("Iteration number: %d\n", n);
        
        // Find solution
        if (CPXmipopt(env, lp)) print_error("Error in find a solution");
        
        // value of objective function
        if (CPXgetobjval (env, lp, &objval)) print_error("Error in CPXgetobjval");
        printf("Objective function value: %lf\n", objval);
        
        // number of variables of the problem
        int ncols = CPXgetnumcols(env, lp);
        
        // final value of variables
        double *xstar = (double *) calloc(ncols, sizeof(double));
        if (CPXgetx(env, lp, xstar, 0, ncols - 1)) print_error("Error in CPXgetx");
        if ( VERBOSE >= 1000) for ( int j = 0; j < ncols; j++ ) printf(" ... qstar[%3d] = %10.2lf \n", j+1, xstar[j]);
        
        // create succ, comp and ncomp of current solution
        build_sol(xstar, inst, succ, comp, &ncomp);
        
        // add SEC constraints
        if(ncomp != 1) add_constraints(inst, env, lp, succ, comp, ncomp, n);
        printf("Number of components: %i\n", ncomp);
        
        t2 = second();
        printf("Time after the iteration: %lf sec.\n\n", t2-t1);
        
        // check if ncomp = 1 in the first phase
        if(first_phase && ncomp == 1) {
            first_phase = 0;
            second_phase = 1;
            ncomp = 9999;
        }
        ++n;
        
        // free solutions
        free(xstar);
    }
    
    // SOLUTION RECAP
    printf("------------ SOLUTION INFO ------------\n");
    printf("z_opt: %f \nIterations: %i \nTSP problem solved in %lf sec.\n", objval, n - 1, t2 - t1);
    printf("---------------------------------------\n");
    
    //print solution
    //print_solution(inst, succ);
    
    free(succ);
    free(comp);
}


// ************************** MODEL DEFINITION ************************** //
//
// obj: compute an optimal solution adding SEC constraints at the end of each run of CPXmipopt
// First version
// Add SEC constraints until ncomp (number of components) is equal to 1
//
// ********************************************************************** //
void loop_method_vers1(instance *inst, CPXENVptr env, CPXLPptr lp) {
    
    int *succ = (int *) calloc(inst->nnodes, sizeof(int));
    int *comp = (int *) calloc(inst->nnodes, sizeof(int));
    // number of components
    int ncomp = 9999;
    // iteration number
    int n = 1;
    // objective function value
    double objval = 0;
    
    printf("Resolve instance \"%s\" with Loop Method\n", inst -> input_file);
    
    // Loop method (Benders)
    while (ncomp != 1) {
        
        printf("Iteration number: %d\n", n);
        
        // Find solution
        if (CPXmipopt(env, lp)) print_error("Error in find a solution");
        
        // value of objective function
        if (CPXgetobjval (env, lp, &objval)) print_error("Error in CPXgetobjval");
        printf("Objective function value: %lf\n", objval);
        
        // number of variables of the problem
        int ncols = CPXgetnumcols(env, lp);
        
        // final value of variables
        double *xstar = (double *) calloc(ncols, sizeof(double));
        if (CPXgetx(env, lp, xstar, 0, ncols - 1)) print_error("Error in CPXgetx");
        if ( VERBOSE >= 1000) for ( int j = 0; j < ncols; j++ ) printf(" ... qstar[%3d] = %10.2lf \n", j+1, xstar[j]);
        
        // create succ, comp and ncomp of current solution
        build_sol(xstar, inst, succ, comp, &ncomp);
        
        // add SEC constraints
        add_constraints(inst, env, lp, succ, comp, ncomp, n);
        printf("Number of components: %i\n\n", ncomp);
        
        if (VERBOSE >= 50) {
            printf("Number of components: %i\n", ncomp);
            for (int i=0; i< inst -> nnodes; i++) {
                printf("Node %2i | succ = %2i | comp = %2i \n", i+1, succ[i]+1, comp[i]);
            }
        }
        ++n;
        
        // free solutions
        free(xstar);
    }
    
    // SOLUTION RECAP
    printf("*** SOLUTION INFO *** \nz_opt: %f \nIterations: %i \n", objval, n - 1);
    
    //print solution
    print_solution(inst, succ);
    
    free(succ);
    free(comp);
}


// add SEC constraints
void add_constraints(instance *inst, CPXENVptr env, CPXLPptr lp, int *succ, int *comp, int ncomp, int n) {
    
    char **cname = (char **) calloc(1, sizeof(char*));
    cname[0] = (char *) calloc(100, sizeof(char));
    char **cname1 = (char **) calloc(1, sizeof(char*));
    cname1[0] = (char *) calloc(100, sizeof(char));
    
    sprintf(cname1[0], "");
    
    int izero = 0;
    
    // max (n * (n - 1) / 2) new constraints
    int *index = (int *) calloc((inst->nnodes * inst->nnodes - 1)/2, sizeof(int));
    double *value = (double *) calloc((inst->nnodes * inst->nnodes - 1)/2, sizeof(double));
    
    int nnz = 0;
    double rhs = 0.0;
    char sense = 'L';
    
    // for each component
    for(int i = 0; i < ncomp; i++) {
        // SEC of "n" iteration and "i+1" component
        sprintf(cname[0], "SEC(%d,%d)", n, i+1);
        nnz = 0;
        rhs = 0.0;
        // Strong subtour elimination for component i
        for(int j = 0; j < inst -> nnodes; j++) {
            if (comp[j] == i + 1) {
                rhs = rhs + 1.0;
                for (int k = j + 1; k < inst -> nnodes; k++ ) {
                    if (comp[k] == i + 1) {
                        index[nnz] = xpos(j, k, inst);
                        value[nnz] = 1.0;
                        ++ nnz;
                    }
                }
            }
        }
        rhs = rhs - 1.0;
        CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, cname1, cname);
    }
    
    // save the model
    if ( VERBOSE >= -100 ) CPXwriteprob(env, lp, "model.lp", NULL);
    
    free(index);
    free(value);
    
    free(cname1[0]);
    free(cname1);
    free(cname[0]);
    free(cname);
}

