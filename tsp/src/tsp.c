//
//  tsp.c
//  cplex
//
//  Created by Matteo Moratello on 13/03/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//
#include "tsp.h"
#include <cut.h>

#define EPS 1e-5

// position
int xpos(int i, int j, instance *inst)
{
    if ( i == j ) print_error(" i == j in xpos" );
    if ( i > j ) return xpos(j,i,inst);
    int pos = i * inst->nnodes + j - (( i + 1 ) * ( i + 2 )) / 2;
    return pos;
}

// distance between two nodes
double dist(int i, int j, instance *inst)
{
    double dist = 0;
    // 2D Euclidean distance
    if ( inst -> euc_2d ) {
        double dx = inst->xcoord[i] - inst->xcoord[j];
        double dy = inst->ycoord[i] - inst->ycoord[j];
        int dis = sqrt(dx*dx+dy*dy) + 0.499999999; // nearest integer
        dist = dis + 0.0;
    }
    // Pseudo-Euclidean distance (only for ATT)
    if ( inst -> att ) {
        double dx = inst->xcoord[i] - inst->xcoord[j];
        double dy = inst->ycoord[i] - inst->ycoord[j];
        double rij = sqrt( (dx*dx + dy*dy) / 10.0 );
        int tij = rij + 0.499999999;               // nearest integer
        int dij = 0;
        if (tij < rij) dij = tij + 1;
        else dij = tij;
        dist = dij + 0.0;
    }
    // Geographical distance
    if ( inst -> geo ) {
        // compute i-th latitude and longitude
        double PI = 3.141592;
        double deg = inst -> xcoord[i] + 0.499999999;
        double min = inst -> xcoord[i] - deg;
        double latitude_i = PI * (deg + 5.0 * min / 3.0 ) / 180.0;
        deg = inst -> ycoord[i] + 0.499999999;
        min = inst -> ycoord[i] - deg;
        double longitude_i = PI * (deg + 5.0 * min / 3.0 ) / 180.0;
        // compute j-th latitude and longitude
        deg = inst -> xcoord[j] + 0.499999999;
        min = inst -> xcoord[j] - deg;
        double latitude_j = PI * (deg + 5.0 * min / 3.0 ) / 180.0;
        deg = inst -> ycoord[j] + 0.499999999;
        min = inst -> ycoord[j] - deg;
        double longitude_j = PI * (deg + 5.0 * min / 3.0 ) / 180.0;
        // compute the distance
        double RRR = 6378.388;
        double q1 = cos( longitude_i - longitude_j );
        double q2 = cos( latitude_i - latitude_j );
        double q3 = cos( latitude_i + latitude_j );
        int dij = (int) ( RRR * acos( 0.5 * ((1.0+q1)*q2 - (1.0-q1)*q3) ) + 1.0);
        dist = dij + 0.0;
    }
    return dist;
}

// optimizer
int TSPopt(instance *inst, double t1) {
    
    // open cplex model
    int error;
    // CPLEX environment
    CPXENVptr env = CPXopenCPLEX(&error);
    // problem object
    CPXLPptr lp = CPXcreateprob(env, &error, "TSP");
    
    //*** CPLEX'S PARAMETERS ***
    CPXsetintparam(env, CPX_PARAM_MIPDISPLAY, 4);
    //CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); // Cplex output on screen
    // set random seed
    CPXsetintparam(env, CPXPARAM_RandomSeed, inst -> randomseed);
    
    //*** CPLEX'S MODEL SETTING ***
    // Standard model
    if(inst -> model_type == 0) build_model_0(inst, env, lp);
    
    // COMPACT MODELS
    // MTZ model
    if(inst -> model_type == 1) build_model_1(inst, env, lp);
    // Flow1 - GG model
    if(inst -> model_type == 2) build_model_2(inst, env, lp);
    // Flow2 - FCG model
    if(inst -> model_type == 3) build_model_3(inst, env, lp);
    // Flow3 - WC model
    if(inst -> model_type == 4) build_model_4(inst, env, lp);
    // T1 - T1 model
    if(inst -> model_type == 5) build_model_5(inst, env, lp);
    // T2 - T2 model
    if(inst -> model_type == 6) build_model_6(inst, env, lp);
    // T3 - T3 model
    if(inst -> model_type == 7) build_model_7(inst, env, lp);
    
    // compact models
    if (inst -> model_type > 0) {
        
        // Find solution
        if (CPXmipopt(env,lp)) print_error("Error in find a solution");
        
        double objval;
        // value of objective function
        if (CPXgetobjval (env, lp, &objval)) print_error("Error in CPXgetobjval");
        printf("Objective function value: %lf\n" , objval);
    
        // number of variables of the problem
        int ncols = CPXgetnumcols(env, lp);
    
        // final value of variables
        double *xstar = (double *) calloc(ncols, sizeof(double));
        if (CPXgetx(env, lp, xstar, 0, ncols-1)) print_error("Error in CPXgetx");
        if ( VERBOSE >= 1000) for ( int j = 0; j < ncols; j++ ) printf(" ... qstar[%3d] = %10.2lf \n", j+1, xstar[j]);
    
        int *succ = (int *) calloc(inst->nnodes, sizeof(int));
        int *comp = (int *) calloc(inst->nnodes, sizeof(int));
        int ncomp = 0;
    
        if(inst -> model_type > 0) build_compact_sol(xstar, inst, succ, comp, &ncomp);
    
        if (VERBOSE >= 50) {
            printf("Number of components: %i\n", ncomp);
            for (int i=0; i< inst -> nnodes; i++) {
                printf("Node %2i | succ = %2i | comp = %2i \n", i+1, succ[i]+1, comp[i]);
            }
        }
    
        // show the solution found
        print_solution(inst, succ);
        
        free(xstar);
        free(succ);
        free(comp);
        
        return 0;
    }
    
    // LOOP METHOD (Benders): loop_method.c
    if (inst -> model_type == 0 && inst -> loop == 1) loop_method(inst, env, lp, t1);
    
    // CALLBACK METHOD: callback.c
    else if(inst -> model_type == 0 && inst -> callback == 1) {
        
        // lazyconstraint callback
        CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);
        
        // First version: pre-generic callbacks
        //if (CPXsetlazyconstraintcallbackfunc(env, mylazycallback, inst)) print_error("Error in lazy constraints callback");
        //if (CPXsetusercutcallbackfunc(env, UserCutCallback, inst)) print_error("Error in user callback");
        // heuristic callback
        //if (CPXsetheuristiccallbackfunc (env, myheuristic, inst)) print_error("Error in heuristic callback");
        
        // Second Version: generic callback
        if (CPXcallbacksetfunc(env, lp, CPX_CALLBACKCONTEXT_CANDIDATE | CPX_CALLBACKCONTEXT_RELAXATION, my_generic_callback, inst)) print_error("Error in generic callback");
        
        int ncores = 1; CPXgetnumcores(env, &ncores);
        CPXsetintparam(env, CPX_PARAM_THREADS, ncores);
        
        inst -> ncols = CPXgetnumcols(env,lp);
        
        // add a starting TSP solution to MIP
        int mcnt = 1;
        int nzcnt = inst -> ncols;
        int beg[1];
        int effort[1];
        beg[0] = 0;
        effort[0] = 4;
        int *varindices = malloc(sizeof(int) * nzcnt);
        double *xstar = (double *) calloc(nzcnt, sizeof(double));
        
        // heuristics
        //NearNeigh(inst, xstar);
        insertion_ch(inst, xstar);
        
        for (int i = 0; i < nzcnt; i++) varindices[i] = i;
        
        if (CPXaddmipstarts(env, lp, mcnt, nzcnt, beg, varindices, xstar, effort, NULL)) print_error("Error in set a mip start");
        
        if (CPXmipopt(env,lp)) print_error("Error in find a solution");
        
        free(xstar);
        free(varindices);
    }
    
    // HARD FIXING
    else if(inst -> model_type == 0 && inst -> hardfixing == 1) hardfixing(inst, env, lp);
    
    // LOCAL BRANCHING
    else if(inst -> model_type == 0 && inst -> localbranching == 1) localbranching(inst, env, lp);
    
    else {
        // Find solution
        if (CPXmipopt(env,lp)) print_error("Error in find a solution");
    }
    
    // *********** Compute and show the final solution *********** //
    
    double objval;
    // value of objective function
    if (CPXgetobjval (env, lp, &objval)) print_error("Error in CPXgetobjval");
    printf("Objective function value: %lf\n" , objval);
    
    // number of variables of the problem
    int ncols = CPXgetnumcols(env, lp);
    
    // final value of variables
    double *xstar = (double *) calloc(ncols, sizeof(double));
    if (CPXgetx(env, lp, xstar, 0, ncols-1)) print_error("Error in CPXgetx");
    
    int *succ = (int *) calloc(inst->nnodes, sizeof(int));
    int *comp = (int *) calloc(inst->nnodes, sizeof(int));
    int ncomp = 0;
    
    build_sol(xstar, inst, succ, comp, &ncomp);
    
    // show the complete solution found (lines + nodes + index node)
    //print_solution(inst, succ);
    // show the solution found (only lines)
    print_solution_light(inst, succ);
    
    free(xstar);
    free(succ);
    free(comp);
    
    // free and close cplex model
    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);
    
    return 0;
}

// build the model
void build_model_0(instance *inst, CPXENVptr env, CPXLPptr lp) {
    
    // type: Binary
    char binary = 'B';
    
    char **cname = (char **) calloc(1, sizeof(char*));
    cname[0] = (char *) calloc(100, sizeof(char));
    
    // add binary var.s x(i,j) for i < j
    for ( int i = 0; i < inst->nnodes; i++ )
    {
        for ( int j = i+1; j < inst->nnodes; j++ )
        {
            sprintf(cname[0], "x(%d,%d)", i+1,j+1);
            double obj = dist(i,j,inst);                // cost == distance
            double lb = 0.0;
            double ub = 1.0;
            if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname) ) print_error(" wrong CPXnewcols on x var.s");
            if ( CPXgetnumcols(env, lp) - 1 != xpos(i, j, inst) ) print_error(" wrong position for x var.s");
        }
    }
    
    int *index = (int *) calloc(inst -> nnodes - 1, sizeof(int));
    double *value = (double *) calloc(inst -> nnodes - 1, sizeof(double));
    int begin[1];
    begin[0] = 0;
    double rhs = 2.0;
    char sense = 'E';                                     // 'E' for equality constraint
    int count = 0;
    
    // degree constraints --> sum{j = 1}^n [x(i,j)] = 2, for each i
    for ( int h = 0; h < inst -> nnodes; h++ )            // degree constraints
    {
        count = 0;
        sprintf(cname[0], "degree(%d)", h+1);
        for ( int i = 0; i < inst->nnodes; i++ ) {
            if ( i == h ) continue;
            index[count] = xpos(i, h, inst);
            value[count++] = 1.0;
        }
        if (CPXaddrows(env, lp, 0, 1, inst -> nnodes - 1, &rhs, &sense, begin, index, value, NULL, cname)) print_error("Error in add row");
    }
    
    // save the model
    if ( VERBOSE >= -100 ) CPXwriteprob(env, lp, "model.lp", NULL);
    
    free(value);
    free(index);
    free(cname[0]);
    free(cname);
    
}

// build succ() and comp() wrt xstar()...
void build_sol(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp) {

// only for debug
#if debug
    int *degree = (int *) calloc(inst->nnodes, sizeof(int));
    for ( int i = 0; i < inst->nnodes; i++ )
    {
        for ( int j = i+1; j < inst->nnodes; j++ )
        {
            int k = xpos(i,j,inst);
            if ( fabs(xstar[k]) > EPS && fabs(xstar[k]-1.0) > EPS ) print_error(" wrong xstar in build_sol()");
            if ( xstar[k] > 0.5 )
            {
                ++degree[i];
                ++degree[j];
            }
        }
    }
    for ( int i = 0; i < inst->nnodes; i++ )
    {
        if ( degree[i] != 2 ) print_error("wrong degree in build_sol()");
    }
    free(degree);
#endif
 
    *ncomp = 0;
    for ( int i = 0; i < inst->nnodes; i++ )
    {
        succ[i] = -1;
        comp[i] = -1;
    }
    
    for ( int start = 0; start < inst->nnodes; start++ )
    {
        if ( comp[start] >= 0 ) continue;  // node "start" was already visited, just skip it
        
        // a new component is found
        (*ncomp)++;
        int i = start;
        while (comp[i] == -1) {
            comp[i] = *ncomp;
            // for each edge
            for ( int j = 0; j < inst->nnodes; j++ ) {
                // check for the first edge of the solution
                if((i != j) && xstar[xpos(i,j,inst)] > 0.5 && comp[j] == -1) {
                    succ[i] = j;
                    i = j;
                    break;
                }
            }
        }
        succ[i] = start;
        // go to the next component...
    }
}
