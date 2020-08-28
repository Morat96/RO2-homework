//
//  compact_models.c
//  cplex
//
//  Created by Matteo Moratello on 13/04/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//

#include "compact_models.h"

#define EPS 1e-5

// **************************** COMPACT MODELS ************************** //
// 1) MTZ
// 2) Flow 1
// 3) Flow 1
// 4) Flow 1
// 5) T1
// 6) T2
// 7) T3
// ********************************************************************** //

/**
 Model 1: compact model MTZ

 @param inst instance of the struct "instance" for TSP problem.
 @param env CPLEX environment.
 @param lp CPLEX LP.
 */
void build_model_1(instance *inst, CPXENVptr env, CPXLPptr lp) {
    
    printf("Resolve instance \"%s\" with Compact Model: MTZ\n\n", inst -> input_file);
    
    // type: Binary
    char binary = 'B';
    // type: Continuous
    char continuous = 'C';
    
    char **cname = (char **) calloc(1, sizeof(char*));
    cname[0] = (char *) calloc(100, sizeof(char));
    
    // ********************************************************************** //
    // ******************************* VARIABLES **************************** //
    // ********************************************************************** //
    
    // add binary var x(i,j)
    for ( int i = 0; i < inst->nnodes; i++ )
    {
        for ( int j = 0; j < inst->nnodes; j++ )
        {
            sprintf(cname[0], "x(%d,%d)", i+1,j+1);
            double obj = dist(i,j,inst); // cost == distance
            double lb = 0.0;
            double ub = ( i == j ) ? 0.0 : 1.0;
            if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname) ) print_error(" wrong CPXnewcols on x var.s");
            if ( CPXgetnumcols(env, lp)-1 != xpos_compact(i, j, inst) ) print_error(" wrong position for x var.s");
        }
    }
    
    // get the last column
    inst -> ustart = CPXgetnumcols(env,lp);
    
    // add continuous var u(i)
    for ( int i = 0; i < inst -> nnodes; i++ )
    {
        sprintf(cname[0], "u(%d)", i+1);
        double obj = 0.0;
        double lb = 0.0;
        double ub = ( i == 0 ) ? 0.0 : inst -> nnodes - 1;
        if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &continuous, cname) ) print_error(" wrong CPXnewcols on u var.s");
        if ( CPXgetnumcols(env,lp)-1 != upos(i, inst) ) print_error(" wrong position for u var.s");
    }
    
    // ********************************************************************** //
    // ***************************** CONSTRAINTS **************************** //
    // ********************************************************************** //
    
    int *cindex = (int *) calloc(inst -> nnodes, sizeof(int));
    double *cvalue = (double *) calloc(inst -> nnodes, sizeof(double));
    int begin[1];
    begin[0] = 0;
    
    // x-constraints (out- and in-degree at nodes)
    for ( int h = 0; h < inst -> nnodes; h++ )  // out-degree
    {
        double rhs = 1.0;
        char sense = 'E';
        sprintf(cname[0], "outdeg(%d)", h+1);
        for ( int i = 0; i < inst->nnodes; i++ )
        {
            cindex[i] = xpos_compact(i, h, inst);
            cvalue[i] = 1.0;
        }
        if (CPXaddrows(env, lp, 0, 1, inst -> nnodes, &rhs, &sense, begin, cindex, cvalue, NULL, cname)) print_error("Error in add row out-deg");
    }
    
    for ( int h = 0; h < inst->nnodes; h++ )   // in-degree
    {
        double rhs = 1.0;
        char sense = 'E';
        sprintf(cname[0], "indeg(%d)", h+1);
        for ( int i = 0; i < inst->nnodes; i++ )
        {
            cindex[i] = xpos_compact(h, i, inst);
            cvalue[i] = 1.0;
        }
        if (CPXaddrows(env, lp, 0, 1, inst -> nnodes, &rhs, &sense, begin, cindex, cvalue, NULL, cname)) print_error("Error in add row in-deg");
    }
    
    free(cindex);
    free(cvalue);
    
    // ********************************************************************** //
    // ************************** LAZY CONSTRAINTS ************************** //
    // ********************************************************************** //
    
    int izero = 0;
    int index[3];
    double value[3];
    
    // add lazy constraints  u(i) - u(j) + M * x(i,j) <= M - 1, for each arc (i,j) not touching node 0
    double big_M = inst -> nnodes;
    double rhs = big_M - 1.0;
    char sense = 'L';
    int nnz = 3;
    
    for ( int i = 1; i < inst->nnodes; i++ ) // excluding node 0
    {
        for ( int j = 1; j < inst->nnodes; j++ ) // excluding node 0
        {
            if ( i == j ) continue;
            sprintf(cname[0], "u-consistency for arc (%d,%d)", i+1, j+1);
            index[0] = upos(i, inst);
            value[0] = 1.0;
            index[1] = upos(j, inst);
            value[1] = - 1.0;
            index[2] = xpos_compact(i, j, inst);
            value[2] = big_M;
            if ( CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero, index, value, cname) ) print_error("wrong CPXlazyconstraints() for u-consistency");
        }
    }
    
    // add lazy constraints x(i,j) + x(j,i) <= 1, for each arc (i,j) with i < j
    rhs = 1.0;
    nnz = 2;
    
    for ( int i = 0; i < inst->nnodes; i++ )
    {
        for ( int j = i+1; j < inst->nnodes; j++ )
        {
            sprintf(cname[0], "SEC on node pair (%d,%d)", i+1, j+1);
            index[0] = xpos_compact(i, j, inst);
            value[0] = 1.0;
            index[1] = xpos_compact(j, i, inst);
            value[1] = 1.0;
            if ( CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero, index, value, cname) ) print_error("wrong CPXlazyconstraints on 2-node SECs");
        }
    }
    
    // save the model
    if ( VERBOSE >= -100 ) CPXwriteprob(env, lp, "model.lp", NULL);
    
    free(cname[0]);
    free(cname);
    
}

/**
 Model 2: compact model FLOW1
 
 @param inst instance of the struct "instance" for TSP problem.
 @param env CPLEX environment.
 @param lp CPLEX LP.
 */
void build_model_2(instance *inst, CPXENVptr env, CPXLPptr lp) {
    
    printf("Resolve instance \"%s\" with Compact Model: Flow1\n\n", inst -> input_file);
    
    // type: Binary
    char binary = 'B';
    // type: Continuous
    char continuous = 'C';
    
    char **cname = (char **) calloc(1, sizeof(char*));
    cname[0] = (char *) calloc(100, sizeof(char));
    
    // ********************************************************************** //
    // ******************************* VARIABLES **************************** //
    // ********************************************************************** //
    
    // add binary var x(i,j) such that 0 <= x(i,j) <= 1
    for ( int i = 0; i < inst -> nnodes; i++ )
    {
        for ( int j = 0; j < inst -> nnodes; j++ )
        {
            sprintf(cname[0], "x(%d,%d)", i+1,j+1);
            double obj = dist(i,j,inst); // cost == distance
            double lb = 0.0;
            double ub = ( i == j ) ? 0.0 : 1.0;
            if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname) ) print_error(" wrong CPXnewcols on x var.s");
            if ( CPXgetnumcols(env,lp)-1 != xpos_compact(i, j, inst) ) print_error(" wrong position for x var.s");
        }
    }
    
    // get the last column
    inst -> ystart = CPXgetnumcols(env, lp);
    
    // add continuous var y(i,j) such that 0 <= y(i,j) <= n
    for ( int i = 0; i < inst -> nnodes; i++ )
    {
        for ( int j = 0; j < inst -> nnodes; j++ )
        {
            sprintf(cname[0], "y(%d,%d)", i+1, j+1);
            double obj = 0.0;
            double lb = 0.0;
            // y(i,i) = 0, for each i
            double ub = ( i == j ) ? 0.0 : inst -> nnodes - 1;
            // y(i,1) = 0, for each i
            if ( j == 0 ) ub = 0.0;
            if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &continuous, cname) ) print_error(" wrong CPXnewcols on y var.s");
            if ( CPXgetnumcols(env,lp)-1 != ypos(i, j, inst) ) print_error(" wrong position for y var.s");
        }
    }
    
    // ********************************************************************** //
    // ***************************** CONSTRAINTS **************************** //
    // ********************************************************************** //
    
    int *cindex = (int *) calloc(inst -> nnodes, sizeof(int));
    double *cvalue = (double *) calloc(inst -> nnodes, sizeof(double));
    int begin[1];
    begin[0] = 0;
    
    // x-constraints (out- and in-degree at nodes)
    for ( int h = 0; h < inst -> nnodes; h++ )  // out-degree
    {
        double rhs = 1.0;
        char sense = 'E';
        sprintf(cname[0], "outdeg(%d)", h+1);
        for ( int i = 0; i < inst->nnodes; i++ )
        {
            cindex[i] = xpos_compact(i, h, inst);
            cvalue[i] = 1.0;
        }
        if (CPXaddrows(env, lp, 0, 1, inst -> nnodes, &rhs, &sense, begin, cindex, cvalue, NULL, cname)) print_error("Error in add row out-deg");
    }
    
    for ( int h = 0; h < inst->nnodes; h++ )   // in-degree
    {
        double rhs = 1.0;
        char sense = 'E';
        sprintf(cname[0], "indeg(%d)", h+1);
        for ( int i = 0; i < inst->nnodes; i++ )
        {
            cindex[i] = xpos_compact(h, i, inst);
            cvalue[i] = 1.0;
        }
        if (CPXaddrows(env, lp, 0, 1, inst -> nnodes, &rhs, &sense, begin, cindex, cvalue, NULL, cname)) print_error("Error in add row in-deg");
    }
    
    free(cindex);
    free(cvalue);
    /*
    // y-constraints --> sum_{j=1}^n [y(j,i)] - sum_{j=1}^n [y(i,j)] = 1, 2 <= i <= n
    // start from y(2,1)
    for ( int h = 1; h < inst->nnodes; h++ )
    {
        int lastrow = CPXgetnumrows(env,lp);
        double rhs = 1.0;
        char sense = 'E';
        sprintf(cname[0], "constr(%d)", h+1);
        if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) print_error(" wrong CPXnewrows [y1]");
        for ( int i = 0; i < inst->nnodes; i++ )
        {
            if (h == i) continue;
            // sum_{n=1}^n [y(j,i)], 2 <= i <= n
            if ( CPXchgcoef(env, lp, lastrow, ypos(i, h, inst), 1.0) ) print_error(" wrong CPXchgcoef [y11]");
            // - sum_{n=1}^n [y(i,j)], 2 <= i <= n
            if ( CPXchgcoef(env, lp, lastrow, ypos(h, i, inst), - 1.0) ) print_error(" wrong CPXchgcoef [y12]");
        }
    }
    */
    // LAZY
    int* ind1 = (int*) calloc(sizeof(int), 2 * inst -> nnodes);
    double* val1 = (double*) calloc(sizeof(double), 2 * inst -> nnodes);
    int izero = 0;
    
    for ( int h = 1; h < inst->nnodes; h++ )
    {
        double rhs = 1.0;
        char sense = 'E';
        sprintf(cname[0], "constr(%d)", h+1);
        int cnt = 0;
        for ( int i = 0; i < inst->nnodes; i++ )
        {
            if (h == i) continue;
            ind1[cnt] = ypos(i, h, inst);
            val1[cnt++] = 1.0;
            ind1[cnt] = ypos(h, i, inst);
            val1[cnt++] = - 1.0;
        }
        if ( CPXaddlazyconstraints(env, lp, 1, cnt, &rhs, &sense, &izero, ind1, val1, cname) ) print_error("wrong CPXlazyconstraints on 2-node SECs");
    }
    
    free(ind1);
    free(val1);
    
    int zero = 0;
    // y-constraints --> sum_{j=1}^n [y(1,j)] = n - 1, j != 1
    int *index = (int *) calloc(inst -> nnodes - 1, sizeof(int));
    double *value = (double *) calloc(inst -> nnodes - 1, sizeof(double));
    
    double rhs = inst -> nnodes - 1;
    char sense = 'E';
    int nnz = inst -> nnodes - 1;
    sprintf(cname[0], "constr1(%d)", 1);
    
    for (int i = 1; i < inst -> nnodes; i++) {
        index[i - 1] = ypos(0, i, inst);
        value[i - 1] = 1.0;
    }
    if ( CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &zero, index, value, cname) ) print_error("wrong CPXlazyconstraints on 2-node SECs");
    //if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, begin, index, value, NULL, cname)) print_error("Error in add row");
    
    // ********************************************************************** //
    // ************************** LAZY CONSTRAINTS ************************** //
    // ********************************************************************** //
    
    int ind[2];
    double val[2];
    //int izero = 0;
    /*
    for (int i=0; i< inst -> nnodes; i++) {
        for(int j=0; j< inst -> nnodes; j++) {
            if (i == j) continue;
            double rhs = 0.0;
            char sense = 'L';
            sprintf(cname[0], "link_lazy(%d,%d)", i+1, j+1);
            ind[0] = xpos_compact(i, j, inst);
            val[0] = - (inst -> nnodes - 1);
            ind[1] = ypos(i, j, inst);
            val[1] = 1.0;
            if ( CPXaddlazyconstraints(env, lp, 1, 2, &rhs, &sense, &izero, ind, val, cname) ) print_error("wrong CPXlazyconstraints on 2-node SECs");
        }
    }
    */
    // y-constraints --> y(i,j) <= (n - 1) * x(i,j), i <= i,j <= n, i != j
    for (int i=0; i< inst -> nnodes; i++) {
        for(int j=0; j< inst->nnodes; j++) {
            if (i == j) continue;
            int lastrow = CPXgetnumrows(env,lp);
            double rhs = 0.0;
            char sense = 'L';
            sprintf(cname[0], "link(%d,%d)", i+1, j+1);
            if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) print_error(" wrong CPXnewrows [xy1]");
            // (n - 1) * x(i,j)
            if ( CPXchgcoef(env, lp, lastrow, xpos_compact(i, j, inst), - (inst -> nnodes - 1)) ) print_error(" wrong CPXchgcoef [xy2]");
            // y(i,j)
            if ( CPXchgcoef(env, lp, lastrow, ypos(i, j, inst), 1.0) ) print_error(" wrong CPXchgcoef [xy3]");
        }
    }
    
    // add lazy constraints x(i,j) + x(j,i) <= 1, for each arc (i,j) with i < j
    rhs = 1.0;
    nnz = 2;
    sense = 'L';
    
    for ( int i = 0; i < inst->nnodes; i++ )
    {
        for ( int j = i+1; j < inst->nnodes; j++ )
        {
            sprintf(cname[0], "SEC on node pair (%d,%d)", i+1, j+1);
            ind[0] = xpos_compact(i, j, inst);
            val[0] = 1.0;
            ind[1] = xpos_compact(j, i, inst);
            val[1] = 1.0;
            if ( CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero, ind, val, cname) ) print_error("wrong CPXlazyconstraints on 2-node SECs");
        }
    }
    
    // save the model
    if ( VERBOSE >= -100 ) CPXwriteprob(env, lp, "model.lp", NULL);
    
    free(index);
    free(value);
    
    free(cname[0]);
    free(cname);
    
}

/**
 Model 3: compact model FLOW2
 
 @param inst instance of the struct "instance" for TSP problem.
 @param env CPLEX environment.
 @param lp CPLEX LP.
 */
void build_model_3(instance *inst, CPXENVptr env, CPXLPptr lp) {
    
    printf("Resolve instance \"%s\" with Compact Model: Flow2\n\n", inst -> input_file);
    
    // type: Binary
    char binary = 'B';
    // type: Continuous
    char continuous = 'C';
    
    char **cname = (char **) calloc(1, sizeof(char*));
    cname[0] = (char *) calloc(100, sizeof(char));
    
    // ********************************************************************** //
    // ******************************* VARIABLES **************************** //
    // ********************************************************************** //
    
    // add binary var x(i,j) such that 0 <= x(i,j) <= 1
    for ( int i = 0; i < inst->nnodes; i++ )
    {
        for ( int j = 0; j < inst->nnodes; j++ )
        {
            sprintf(cname[0], "x(%d,%d)", i+1,j+1);
            double obj = dist(i,j,inst); // cost == distance
            double lb = 0.0;
            double ub = ( i == j ) ? 0.0 : 1.0;
            if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname) ) print_error(" wrong CPXnewcols on x var.s");
            if ( CPXgetnumcols(env,lp)-1 != xpos_compact(i, j, inst) ) print_error(" wrong position for x var.s");
        }
    }
    
    // get the last column
    inst -> ystart = CPXgetnumcols(env,lp);
    
    // add continuous var y(i,j) such that 0 <= y(i,j) <= n - 1
    for ( int i = 0; i < inst -> nnodes; i++ )
    {
        for ( int j = 0; j < inst -> nnodes; j++ )
        {
            sprintf(cname[0], "y(%d,%d)", i+1, j+1);
            double obj = 0.0;
            double lb = 0.0;
            double ub = ( i == j ) ? 0.0 : inst -> nnodes - 1;
            if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &continuous, cname) ) print_error(" wrong CPXnewcols on y var.s");
            if ( CPXgetnumcols(env,lp) - 1 != ypos(i, j, inst) ) print_error(" wrong position for y var.s");
        }
    }
    
    // get the last column
    inst -> zstart = CPXgetnumcols(env,lp);
    
    // add continuous var y(i,j) such that 0 <= z(i,j) <= n - 1
    for ( int i = 0; i < inst -> nnodes; i++ )
    {
        for ( int j = 0; j < inst -> nnodes; j++ )
        {
            sprintf(cname[0], "z(%d,%d)", i+1, j+1);
            double obj = 0.0;
            double lb = 0.0;
            double ub = ( i == j ) ? 0.0 : inst -> nnodes - 1;
            if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &continuous, cname) ) print_error(" wrong CPXnewcols on z var.s");
            if ( CPXgetnumcols(env,lp)-1 != zpos(i, j, inst) ) print_error(" wrong position for z var.s");
        }
    }
    
    // ********************************************************************** //
    // ***************************** CONSTRAINTS **************************** //
    // ********************************************************************** //
    
    int *cindex = (int *) calloc(inst -> nnodes, sizeof(int));
    double *cvalue = (double *) calloc(inst -> nnodes, sizeof(double));
    int begin[1];
    begin[0] = 0;
    
    // x-constraints (out- and in-degree at nodes)
    for ( int h = 0; h < inst -> nnodes; h++ )  // out-degree
    {
        double rhs = 1.0;
        char sense = 'E';
        sprintf(cname[0], "outdeg(%d)", h+1);
        for ( int i = 0; i < inst->nnodes; i++ )
        {
            cindex[i] = xpos_compact(i, h, inst);
            cvalue[i] = 1.0;
        }
        if (CPXaddrows(env, lp, 0, 1, inst -> nnodes, &rhs, &sense, begin, cindex, cvalue, NULL, cname)) print_error("Error in add row out-deg");
    }
    
    for ( int h = 0; h < inst->nnodes; h++ )   // in-degree
    {
        double rhs = 1.0;
        char sense = 'E';
        sprintf(cname[0], "indeg(%d)", h+1);
        for ( int i = 0; i < inst->nnodes; i++ )
        {
            cindex[i] = xpos_compact(h, i, inst);
            cvalue[i] = 1.0;
        }
        if (CPXaddrows(env, lp, 0, 1, inst -> nnodes, &rhs, &sense, begin, cindex, cvalue, NULL, cname)) print_error("Error in add row in-deg");
    }
    
    free(cindex);
    free(cvalue);
    
    // ********************************************************************** //
    // ************************** LAZY CONSTRAINTS ************************** //
    // ********************************************************************** //
    
    int izero = 0;
    int *index = (int *) calloc(2 * (inst -> nnodes - 1), sizeof(int));
    double *value = (double *) calloc(2 * (inst -> nnodes - 1), sizeof(double));
    
    // first y-constraints --> sum_{j=1,j != 1}^n [y(1,j) - y(j,1)] = n - 1,
    double rhs = inst -> nnodes - 1;
    char sense = 'E';
    int nnz = 2 * (inst -> nnodes - 1);
    
    sprintf(cname[0], "constr_y_1(%d)", 1);
    for (int i = 1; i < inst -> nnodes; i++) {
        index[i - 1] = ypos(0, i, inst);
        value[i - 1] = 1.0;
        index[inst -> nnodes + i - 2] = ypos(0, i * (inst -> nnodes), inst);
        value[inst -> nnodes + i - 2] = 1.0;
    }
    //if ( CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero, index, value, cname) ) print_error("wrong CPXlazyconstraints on 2-node SECs");
    if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, begin, index, value, NULL, cname)) print_error("Error in add row");
    
    // first y-constraints --> sum_{j=1,j != 1}^n [y(1,j) - y(j,1)] = n - 1,
    sprintf(cname[0], "constr_z_1(%d)", 1);
    for (int i = 1; i < inst -> nnodes; i++) {
        index[i - 1] = zpos(0, i, inst);
        value[i - 1] = 1.0;
        index[inst -> nnodes + i - 2] = zpos(0, i * (inst -> nnodes), inst);
        value[inst -> nnodes + i - 2] = 1.0;
    }
    //if ( CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero, index, value, cname) ) print_error("wrong CPXlazyconstraints on 2-node SECs");
    if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, begin, index, value, NULL, cname)) print_error("Error in add row");
    
    free(index);
    free(value);
    
    // ********************************************************************** //
    
    int* ind1 = (int*) calloc(sizeof(int), 2 * inst -> nnodes);
    double* val1 = (double*) calloc(sizeof(double), 2 * inst -> nnodes);
    
    // second y-constraints --> sum_{j=1}^n [y(i,j) - y(j,i)] = -1, for each i\{1}, i != j
    for ( int h = 1; h < inst->nnodes; h++ )
    {
        int lastrow = CPXgetnumrows(env,lp);
        double rhs = - 1.0;
        char sense = 'E';
        sprintf(cname[0], "constr_y(%d)", h+1);
        int cnt = 0;
        if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) print_error(" wrong CPXnewrows [y2]");
        for ( int i = 0; i < inst->nnodes; i++ )
        {
            if (h == i) continue;
            // y(j,i)
            
            if ( CPXchgcoef(env, lp, lastrow, ypos(h, i, inst), 1.0) ) print_error(" wrong CPXchgcoef [y21]");
            // - y(j,i)
            if ( CPXchgcoef(env, lp, lastrow, ypos(i, h, inst), - 1.0)) print_error(" wrong CPXchgcoef [y22]");
            
            // LAZY
            /*
            ind1[cnt] = ypos(h, i, inst);
            val1[cnt++] = 1.0;
            ind1[cnt] = ypos(i, h, inst);
            val1[cnt++] = - 1.0;*/
        }
        //if ( CPXaddlazyconstraints(env, lp, 1, cnt, &rhs, &sense, &izero, ind1, val1, cname) ) print_error("wrong CPXlazyconstraints on 2-node SECs");
    }
    
    // second z-constraints --> sum_{j=1}^n [z(i,j) - z(j,i)] = 1, for each i\{1}, i != j
    for ( int h = 1; h < inst->nnodes; h++ )
    {
        int lastrow = CPXgetnumrows(env,lp);
        double rhs = 1.0;
        char sense = 'E';
        sprintf(cname[0], "constr_z(%d)", h+1);
        int cnt = 0;
        //if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) print_error(" wrong CPXnewrows [z2]");
        for ( int i = 0; i < inst->nnodes; i++ )
        {
            
            if (h == i) continue;
            // z(i,j)
            /*
            if ( CPXchgcoef(env, lp, lastrow, zpos(h, i, inst), 1.0) ) print_error(" wrong CPXchgcoef [z21]");
            // - z(j,i)
            if ( CPXchgcoef(env, lp, lastrow, zpos(i, h, inst), - 1.0)) print_error(" wrong CPXchgcoef [z22]");
            */
            
            ind1[cnt] = zpos(h, i, inst);
            val1[cnt++] = 1.0;
            ind1[cnt] = zpos(i, h, inst);
            val1[cnt++] = - 1.0;
        }
        if ( CPXaddlazyconstraints(env, lp, 1, cnt, &rhs, &sense, &izero, ind1, val1, cname) ) print_error("wrong CPXlazyconstraints on 2-node SECs");
    }
    

    
    // y and z links --> sum_{j=1}^n [y(i,j) + z(i,j)] = n - 1, for each i
    for ( int h = 0; h < inst->nnodes; h++ )
    {
        int lastrow = CPXgetnumrows(env,lp);
        double rhs = inst -> nnodes - 1;
        char sense = 'E';
        sprintf(cname[0], "link_y_z(%d)", h+1);
        int cnt = 0;
        if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) print_error(" wrong CPXnewrows [yz1]");
        for ( int i = 0; i < inst->nnodes; i++ )
        {
            
            // y(i,j)
            if ( CPXchgcoef(env, lp, lastrow, ypos(h, i, inst), 1.0) ) print_error(" wrong CPXchgcoef [yz1]");
            // z(i,j)
            if ( CPXchgcoef(env, lp, lastrow, zpos(h, i, inst), 1.0)) print_error(" wrong CPXchgcoef [yz2]");
            /*
            ind1[cnt] = ypos(h, i, inst);
            val1[cnt++] = 1.0;
            ind1[cnt] = zpos(h, i, inst);
            val1[cnt++] = 1.0;*/
        }
        //if ( CPXaddlazyconstraints(env, lp, 1, cnt, &rhs, &sense, &izero, ind1, val1, cname) ) print_error("wrong CPXlazyconstraints on 2-node SECs");
    }
    
    free(ind1);
    free(val1);
    
    int ind[3];
    double val[3];
    
    // x, y and z links --> y(i,j) + z(i,j) = (n - 1) * x(i,j), for each i,j
    for ( int h = 0; h < inst->nnodes; h++ )
    {
        for ( int i = 0; i < inst->nnodes; i++ )
        {
            
            int lastrow = CPXgetnumrows(env,lp);
            double rhs = 0.0;
            char sense = 'E';
            sprintf(cname[0], "link_x_y_z(%d,%d)", h+1, i+1);
            if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) print_error(" wrong CPXnewrows [xyz1]");
            // y(i,j)
            if ( CPXchgcoef(env, lp, lastrow, ypos(h, i, inst), 1.0) ) print_error(" wrong CPXchgcoef [xyz1]");
            // z(i,j)
            if ( CPXchgcoef(env, lp, lastrow, zpos(h, i, inst), 1.0)) print_error(" wrong CPXchgcoef [xyz2]");
            // (n - 1) * x(i,j)
            if ( CPXchgcoef(env, lp, lastrow, xpos_compact(h, i, inst), - (inst -> nnodes -1))) print_error(" wrong CPXchgcoef [xyz3]");
            /*ind[0] = ypos(h, i, inst);
            val[0] = 1.0;
            ind[1] = zpos(h, i, inst);
            val[1] = 1.0;
            ind[2] = xpos_compact(h, i, inst);
            val[2] = - (inst -> nnodes - 1);*/
            
            //if ( CPXaddlazyconstraints(env, lp, 1, 3, &rhs, &sense, &izero, ind, val, cname) ) print_error("wrong CPXlazyconstraints on 2-node SECs");
            
            
        }
    }
    
    // save the model
    if ( VERBOSE >= -100 ) CPXwriteprob(env, lp, "model.lp", NULL);
    
    free(cname[0]);
    free(cname);
    
}

/**
 Model 4: compact model FLOW3
 
 @param inst instance of the struct "instance" for TSP problem.
 @param env CPLEX environment.
 @param lp CPLEX LP.
 */
void build_model_4(instance *inst, CPXENVptr env, CPXLPptr lp) {
    
    printf("Resolve instance \"%s\" with Compact Model: Flow3\n\n", inst -> input_file);
    
    // type: Binary
    char binary = 'B';
    // type: Continuous
    char continuous = 'C';
    
    char **cname = (char **) calloc(1, sizeof(char*));
    cname[0] = (char *) calloc(100, sizeof(char));
    
    // ********************************************************************** //
    // ******************************* VARIABLES **************************** //
    // ********************************************************************** //
    
    // add binary var x(i,j) such that 0 <= x(i,j) <= 1
    for ( int i = 0; i < inst->nnodes; i++ )
    {
        for ( int j = 0; j < inst->nnodes; j++ )
        {
            sprintf(cname[0], "x(%d,%d)", i+1,j+1);
            double obj = dist(i,j,inst); // cost == distance
            double lb = 0.0;
            double ub = ( i == j ) ? 0.0 : 1.0;
            if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname) ) print_error(" wrong CPXnewcols on x var.s");
            if ( CPXgetnumcols(env,lp)-1 != xpos_compact(i, j, inst) ) print_error(" wrong position for x var.s");
        }
    }
    
    // get the last column
    inst -> ystart = CPXgetnumcols(env, lp);
    
    // add continuous var y(i,j) such that 0 <= y(i,j) <= n
    for (int k = 1; k < inst -> nnodes; k++) {
        for ( int i = 0; i < inst -> nnodes; i++ ) {
            for ( int j = 0; j < inst -> nnodes; j++ ) {
                sprintf(cname[0], "y(%d,%d,%d)",k+1, i+1, j+1);
                double obj = 0.0;
                double lb = 0.0;
                double ub = inst -> nnodes - 1;
                if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &continuous, cname) ) print_error(" wrong CPXnewcols on y var.s");
                if ( CPXgetnumcols(env,lp) - 1 != y2pos(i, j, k, inst) ) print_error(" wrong position for y var.s");
            }
        }
    }
    
    // ********************************************************************** //
    // ***************************** CONSTRAINTS **************************** //
    // ********************************************************************** //
    
    int *cindex = (int *) calloc(inst -> nnodes, sizeof(int));
    double *cvalue = (double *) calloc(inst -> nnodes, sizeof(double));
    int begin[1];
    begin[0] = 0;
    
    // x-constraints (out- and in-degree at nodes)
    for ( int h = 0; h < inst -> nnodes; h++ )  // out-degree
    {
        double rhs = 1.0;
        char sense = 'E';
        sprintf(cname[0], "outdeg(%d)", h+1);
        for ( int i = 0; i < inst->nnodes; i++ )
        {
            cindex[i] = xpos_compact(i, h, inst);
            cvalue[i] = 1.0;
        }
        if (CPXaddrows(env, lp, 0, 1, inst -> nnodes, &rhs, &sense, begin, cindex, cvalue, NULL, cname)) print_error("Error in add row out-deg");
    }
    
    for ( int h = 0; h < inst->nnodes; h++ )   // in-degree
    {
        double rhs = 1.0;
        char sense = 'E';
        sprintf(cname[0], "indeg(%d)", h+1);
        for ( int i = 0; i < inst->nnodes; i++ )
        {
            cindex[i] = xpos_compact(h, i, inst);
            cvalue[i] = 1.0;
        }
        if (CPXaddrows(env, lp, 0, 1, inst -> nnodes, &rhs, &sense, begin, cindex, cvalue, NULL, cname)) print_error("Error in add row in-deg");
    }
    
    free(cindex);
    free(cvalue);
    
    // y(i,j,k) <= x(i,j) for each i,j,k, k != 1
    for (int k = 1; k < inst -> nnodes; k++) {
        for ( int i = 0; i < inst -> nnodes; i++ ) {
            for ( int j = 0; j < inst -> nnodes; j++ ) {
                int lastrow = CPXgetnumrows(env,lp);
                double rhs = 0.0;
                char sense = 'L';
                sprintf(cname[0], "constr(%d,%d,%d)", k+1, i+1, j+1);
                if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) print_error(" wrong CPXnewrows [y1]");
                if ( CPXchgcoef(env, lp, lastrow, xpos_compact(i, j, inst), - 1.0) ) print_error(" wrong CPXchgcoef [y11]");
                if ( CPXchgcoef(env, lp, lastrow, y2pos(i, j, k, inst), 1.0) ) print_error(" wrong CPXchgcoef [y12]");
                
            }
        }
    }
    
    // sum {i=1}^n [y(1,i,k)] = 1, for each k - {1}
    for (int k = 1; k < inst -> nnodes; k++) {
        int lastrow = CPXgetnumrows(env,lp);
        double rhs = 1.0;
        char sense = 'E';
        sprintf(cname[0], "constr(%d)", k+1);
        if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) print_error(" wrong CPXnewrows [y1]");
        for ( int i = 0; i < inst -> nnodes; i++ ) {
            if ( CPXchgcoef(env, lp, lastrow, y2pos(0, i, k, inst), 1.0) ) print_error(" wrong CPXchgcoef [y11]");
        }
    }
    
    // sum {i=1}^n [y(i,1,k)] = 0, for each k - {1}
    for (int k = 1; k < inst -> nnodes; k++) {
        int lastrow = CPXgetnumrows(env,lp);
        double rhs = 0.0;
        char sense = 'E';
        sprintf(cname[0], "constr2(%d)", k+1);
        if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) print_error(" wrong CPXnewrows [y1]");
        for ( int i = 0; i < inst -> nnodes; i++ ) {
            if ( CPXchgcoef(env, lp, lastrow, y2pos(i, 0, k, inst), 1.0) ) print_error(" wrong CPXchgcoef [y11]");
        }
    }
    
    // sum {i=1}^n [y(i,k,k)] = 1, for each k - {1}
    for (int k = 1; k < inst -> nnodes; k++) {
        int lastrow = CPXgetnumrows(env,lp);
        double rhs = 1.0;
        char sense = 'E';
        sprintf(cname[0], "constr3(%d)", k+1);
        if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) print_error(" wrong CPXnewrows [y1]");
        for ( int i = 0; i < inst -> nnodes; i++ ) {
            if ( CPXchgcoef(env, lp, lastrow, y2pos(i, k, k, inst), 1.0) ) print_error(" wrong CPXchgcoef [y11]");
        }
    }
    
    // sum {j=1}^n [y(k,j,k)] = 0, for each k - {1}
    for (int k = 1; k < inst -> nnodes; k++) {
        int lastrow = CPXgetnumrows(env,lp);
        double rhs = 0.0;
        char sense = 'E';
        sprintf(cname[0], "constr4(%d)", k+1);
        if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) print_error(" wrong CPXnewrows [y1]");
        for ( int i = 0; i < inst -> nnodes; i++ ) {
            if ( CPXchgcoef(env, lp, lastrow, y2pos(k, i, k, inst), 1.0) ) print_error(" wrong CPXchgcoef [y11]");
        }
    }
    
    // sum {i=1}^n [y(i,j,k)] - sum {i=1}^n [y(j,i,k)] = 0 for each j,k - {1}, j |= k
    for (int k = 1; k < inst -> nnodes; k++) {
        for (int j = 1; j < inst -> nnodes; j++) {
            if ( k == j ) continue;
            int lastrow = CPXgetnumrows(env,lp);
            double rhs = 0.0;
            char sense = 'E';
            sprintf(cname[0], "constr(%d,%d)", k+1, j+1);
            if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) print_error(" wrong CPXnewrows [y1]");
            for ( int i = 0; i < inst -> nnodes; i++ ) {
                if ( CPXchgcoef(env, lp, lastrow, y2pos(i, j, k, inst), 1.0) ) print_error(" wrong CPXchgcoef [y11]");
                if ( CPXchgcoef(env, lp, lastrow, y2pos(j, i, k, inst), - 1.0) ) print_error(" wrong CPXchgcoef [y11]");
            }
        }
    }
    
    // save the model
    if ( VERBOSE >= -100 ) CPXwriteprob(env, lp, "model.lp", NULL);
    
    free(cname[0]);
    free(cname);
    
}

/**
 Model 5: compact model T1
 
 @param inst instance of the struct "instance" for TSP problem.
 @param env CPLEX environment.
 @param lp CPLEX LP.
 */
void build_model_5(instance *inst, CPXENVptr env, CPXLPptr lp) {
    
    printf("Resolve instance \"%s\" with Compact Model: T1\n\n", inst -> input_file);
    
    // type: Binary
    char binary = 'B';
    // type: Continuous
    char continuous = 'C';
    
    char **cname = (char **) calloc(1, sizeof(char*));
    cname[0] = (char *) calloc(100, sizeof(char));
    
    // ********************************************************************** //
    // ******************************* VARIABLES **************************** //
    // ********************************************************************** //
    
    // add binary var x(i,j) such that 0 <= x(i,j) <= 1
    for ( int i = 0; i < inst->nnodes; i++ )
    {
        for ( int j = 0; j < inst->nnodes; j++ )
        {
            sprintf(cname[0], "x(%d,%d)", i+1,j+1);
            double obj = dist(i,j,inst); // cost == distance
            double lb = 0.0;
            double ub = ( i == j ) ? 0.0 : 1.0;
            if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname) ) print_error(" wrong CPXnewcols on x var.s");
            if ( CPXgetnumcols(env,lp)-1 != xpos_compact(i, j, inst) ) print_error(" wrong position for x var.s");
        }
    }
    
    // get the last column
    inst -> ystart = CPXgetnumcols(env, lp);
    
    // add continuous var y(i,j) such that 0 <= y(i,j) <= n
    for (int t = 0; t < inst -> nnodes; t++) {
        for ( int i = 0; i < inst -> nnodes; i++ ) {
            for ( int j = 0; j < inst -> nnodes; j++ ) {
                sprintf(cname[0], "y(%d,%d,%d)",t+1, i+1, j+1);
                double obj = 0.0;
                double lb = 0.0;
                double ub = inst -> nnodes - 1;
                if (j == 0 && t != (inst -> nnodes - 1)) ub = 0;
                if (i == 0 && t != 0) ub = 0;
                if (t == 0 && i != 0) ub = 0;
                if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &continuous, cname) ) print_error(" wrong CPXnewcols on y var.s");
                if ( CPXgetnumcols(env,lp) - 1 != y3pos(i, j, t, inst) ) print_error(" wrong position for y var.s");
            }
        }
    }
    
    // ********************************************************************** //
    // ***************************** CONSTRAINTS **************************** //
    // ********************************************************************** //
    
    int *cindex = (int *) calloc(inst -> nnodes, sizeof(int));
    double *cvalue = (double *) calloc(inst -> nnodes, sizeof(double));
    int begin[1];
    begin[0] = 0;
    
    // x-constraints (out- and in-degree at nodes)
    for ( int h = 0; h < inst -> nnodes; h++ )  // out-degree
    {
        double rhs = 1.0;
        char sense = 'E';
        sprintf(cname[0], "outdeg(%d)", h+1);
        for ( int i = 0; i < inst->nnodes; i++ )
        {
            cindex[i] = xpos_compact(i, h, inst);
            cvalue[i] = 1.0;
        }
        if (CPXaddrows(env, lp, 0, 1, inst -> nnodes, &rhs, &sense, begin, cindex, cvalue, NULL, cname)) print_error("Error in add row out-deg");
    }
    
    for ( int h = 0; h < inst->nnodes; h++ )   // in-degree
    {
        double rhs = 1.0;
        char sense = 'E';
        sprintf(cname[0], "indeg(%d)", h+1);
        for ( int i = 0; i < inst->nnodes; i++ )
        {
            cindex[i] = xpos_compact(h, i, inst);
            cvalue[i] = 1.0;
        }
        if (CPXaddrows(env, lp, 0, 1, inst -> nnodes, &rhs, &sense, begin, cindex, cvalue, NULL, cname)) print_error("Error in add row in-deg");
    }
    
    free(cindex);
    free(cvalue);
    
    // sum {i,j,t}^n [y(i,j,t)] = n
    int lastrow = CPXgetnumrows(env,lp);
    double rhs = inst -> nnodes;
    char sense = 'E';
    sprintf(cname[0], "constr(%d)", 1);
    if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) print_error(" wrong CPXnewrows [y1]");
    for (int t = 0; t < inst -> nnodes; t++) {
        for ( int i = 0; i < inst -> nnodes; i++ ) {
            for ( int j = 0; j < inst -> nnodes; j++ ) {
                if ( CPXchgcoef(env, lp, lastrow, y3pos(i, j, t, inst), 1.0) ) print_error(" wrong CPXchgcoef [y11]");
            }
        }
    }
    
    // sum {j,t >=2}^n [t * y(i,j,t)] - sum {k,t}^n [t * y(k,i,t)] = 1
    for (int i = 1; i < inst -> nnodes; i++) {
        int lastrow = CPXgetnumrows(env,lp);
        double rhs = 1.0;
        char sense = 'E';
        sprintf(cname[0], "constr1(%d)", i+1);
        if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) print_error(" wrong CPXnewrows [y1]");
        for ( int j = 0; j < inst -> nnodes; j++ ) {
            for ( int k = 1; k < inst -> nnodes; k++ ) {
                if ( CPXchgcoef(env, lp, lastrow, y3pos(i, j, k, inst), k) ) print_error(" wrong CPXchgcoef [y11]");
                if ( CPXchgcoef(env, lp, lastrow, y3pos(j, i, k, inst), - k) ) print_error(" wrong CPXchgcoef [y11]");
            }
        }
        for ( int j = 1; j < inst -> nnodes; j++ ) {
            if ( CPXchgcoef(env, lp, lastrow, y3pos(j, i, 0, inst), - 1) ) print_error(" wrong CPXchgcoef [y11]");
        }
    }
    
    // x(i,j) - sum {t}^n [y(i,j,t)] = 0
    for (int i = 0; i < inst -> nnodes; i++) {
        for (int j = 0; j < inst -> nnodes; j++) {
            if( i == j ) continue;
            int lastrow = CPXgetnumrows(env,lp);
            double rhs = 0.0;
            char sense = 'E';
            sprintf(cname[0], "constr(%d,%d)", i+1, j+1);
            if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) print_error(" wrong CPXnewrows [y1]");
            if ( CPXchgcoef(env, lp, lastrow, xpos_compact(i, j, inst), 1.0) ) print_error(" wrong CPXchgcoef [y11]");
            for ( int t = 0; t < inst -> nnodes; t++ ) {
                if ( CPXchgcoef(env, lp, lastrow, y3pos(i, j, t, inst), - 1.0) ) print_error(" wrong CPXchgcoef [y11]");
            }
        }
    }
    
    // ********************************************************************** //
    // ************************** LAZY CONSTRAINTS ************************** //
    // ********************************************************************** //
    
    int izero = 0;
    int index[2];
    double value[2];
    rhs = 0.0;
    sense = 'E';
    
    // y(i,j,1) - y(j,i,n) = 0, for each i != 1, j
    for (int i = 1; i < inst -> nnodes; i++) {
        for (int j = 0; j < inst -> nnodes; j++) {
            sprintf(cname[0], "constr1(%d,%d)", i+1, j+1);
            index[0] = y3pos(i, j, 0, inst);
            value[0] = 1.0;
            index[1] = y3pos(j, i, inst -> nnodes - 1, inst);
            value[1] = - 1.0;
            if ( CPXaddlazyconstraints(env, lp, 1, 2, &rhs, &sense, &izero, index, value, cname) ) print_error("wrong CPXlazyconstraints on 2-node SECs");
        }
    }
    
    // without lazy constraint
    /*
     // sum {i=1}^n [y(i,1,k)] = 0, for each k - {1}
     for (int i = 1; i < inst -> nnodes; i++) {
     for (int j = 0; j < inst -> nnodes; j++) {
     int lastrow = CPXgetnumrows(env,lp);
     double rhs = 0.0;
     char sense = 'E';
     sprintf(cname[0], "constr3(%d,%d)", i+1, j+1);
     if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) print_error(" wrong CPXnewrows [y1]");
     if ( CPXchgcoef(env, lp, lastrow, y3pos(i, j, 0, inst), 1.0) ) print_error(" wrong CPXchgcoef [y11]");
     if ( CPXchgcoef(env, lp, lastrow, y3pos(j, i, inst -> nnodes - 1, inst), - 1.0) ) print_error(" wrong CPXchgcoef [y11]");
     }
     }*/
    
    // save the model
    if ( VERBOSE >= -100 ) CPXwriteprob(env, lp, "model.lp", NULL);
    
    free(cname[0]);
    free(cname);
    
}

/**
 Model 6: compact model T2
 
 @param inst instance of the struct "instance" for TSP problem.
 @param env CPLEX environment.
 @param lp CPLEX LP.
 */
void build_model_6(instance *inst, CPXENVptr env, CPXLPptr lp) {
    
    printf("Resolve instance \"%s\" with Compact Model: T2\n\n", inst -> input_file);
    
    // type: Binary
    char binary = 'B';
    // type: Continuous
    char continuous = 'C';
    
    char **cname = (char **) calloc(1, sizeof(char*));
    cname[0] = (char *) calloc(100, sizeof(char));
    
    // ********************************************************************** //
    // ******************************* VARIABLES **************************** //
    // ********************************************************************** //
    
    // add binary var x(i,j) such that 0 <= x(i,j) <= 1
    for ( int i = 0; i < inst->nnodes; i++ )
    {
        for ( int j = 0; j < inst->nnodes; j++ )
        {
            sprintf(cname[0], "x(%d,%d)", i+1,j+1);
            double obj = dist(i,j,inst); // cost == distance
            double lb = 0.0;
            double ub = ( i == j ) ? 0.0 : 1.0;
            if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname) ) print_error(" wrong CPXnewcols on x var.s");
            if ( CPXgetnumcols(env,lp)-1 != xpos_compact(i, j, inst) ) print_error(" wrong position for x var.s");
        }
    }
    
    // get the last column
    inst -> ystart = CPXgetnumcols(env, lp);
    
    // add continuous var y(i,j) such that 0 <= y(i,j) <= n
    for (int t = 0; t < inst -> nnodes; t++) {
        for ( int i = 0; i < inst -> nnodes; i++ ) {
            for ( int j = 0; j < inst -> nnodes; j++ ) {
                sprintf(cname[0], "y(%d,%d,%d)",t+1, i+1, j+1);
                double obj = 0.0;
                double lb = 0.0;
                double ub = inst -> nnodes - 1;
                if (j == 0 && t != (inst -> nnodes - 1)) ub = 0;
                if (i == 0 && t != 0) ub = 0;
                if (t == 0 && i != 0) ub = 0;
                if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &continuous, cname) ) print_error(" wrong CPXnewcols on y var.s");
                if ( CPXgetnumcols(env,lp) - 1 != y3pos(i, j, t, inst) ) print_error(" wrong position for y var.s");
            }
        }
    }
    
    // ********************************************************************** //
    // ***************************** CONSTRAINTS **************************** //
    // ********************************************************************** //
    
    int *cindex = (int *) calloc(inst -> nnodes, sizeof(int));
    double *cvalue = (double *) calloc(inst -> nnodes, sizeof(double));
    int begin[1];
    begin[0] = 0;
    
    // x-constraints (out- and in-degree at nodes)
    for ( int h = 0; h < inst -> nnodes; h++ )  // out-degree
    {
        double rhs = 1.0;
        char sense = 'E';
        sprintf(cname[0], "outdeg(%d)", h+1);
        for ( int i = 0; i < inst->nnodes; i++ )
        {
            cindex[i] = xpos_compact(i, h, inst);
            cvalue[i] = 1.0;
        }
        if (CPXaddrows(env, lp, 0, 1, inst -> nnodes, &rhs, &sense, begin, cindex, cvalue, NULL, cname)) print_error("Error in add row out-deg");
    }
    
    for ( int h = 0; h < inst->nnodes; h++ )   // in-degree
    {
        double rhs = 1.0;
        char sense = 'E';
        sprintf(cname[0], "indeg(%d)", h+1);
        for ( int i = 0; i < inst->nnodes; i++ )
        {
            cindex[i] = xpos_compact(h, i, inst);
            cvalue[i] = 1.0;
        }
        if (CPXaddrows(env, lp, 0, 1, inst -> nnodes, &rhs, &sense, begin, cindex, cvalue, NULL, cname)) print_error("Error in add row in-deg");
    }
    
    free(cindex);
    free(cvalue);
    
    // x(i,j) - sum {t}^n [y(i,j,t)] = 0
    for (int i = 0; i < inst -> nnodes; i++) {
        for (int j = 0; j < inst -> nnodes; j++) {
            if( i == j ) continue;
            int lastrow = CPXgetnumrows(env,lp);
            double rhs = 0.0;
            char sense = 'E';
            sprintf(cname[0], "constr(%d,%d)", i+1, j+1);
            if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) print_error(" wrong CPXnewrows [y1]");
            if ( CPXchgcoef(env, lp, lastrow, xpos_compact(i, j, inst), 1.0) ) print_error(" wrong CPXchgcoef [y11]");
            for ( int t = 0; t < inst -> nnodes; t++ ) {
                if ( CPXchgcoef(env, lp, lastrow, y3pos(i, j, t, inst), - 1.0) ) print_error(" wrong CPXchgcoef [y11]");
            }
        }
    }
    
    // sum {i != j,t}^n [y(i,j,t)] = 1, for each j
    for (int j = 0; j < inst -> nnodes; j++) {
        int lastrow = CPXgetnumrows(env,lp);
        double rhs = 1.0;
        char sense = 'E';
        sprintf(cname[0], "constr(%d)", j+1);
        if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) print_error(" wrong CPXnewrows [y1]");
        for (int t = 0; t < inst -> nnodes; t++) {
            for (int i = 0; i < inst -> nnodes; i++) {
                if (i == j) continue;
                if ( CPXchgcoef(env, lp, lastrow, y3pos(i, j, t, inst), 1.0) ) print_error(" wrong CPXchgcoef [y11]");
            }
        }
    }
    
    // sum {j != i,t}^n [y(i,j,t)] = 1, for each i
    for (int i = 0; i < inst -> nnodes; i++) {
        int lastrow = CPXgetnumrows(env,lp);
        double rhs = 1.0;
        char sense = 'E';
        sprintf(cname[0], "constr1(%d)", i+1);
        if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) print_error(" wrong CPXnewrows [y1]");
        for (int t = 0; t < inst -> nnodes; t++) {
            for (int j = 0; j < inst -> nnodes; j++) {
                if (j == i) continue;
                if ( CPXchgcoef(env, lp, lastrow, y3pos(i, j, t, inst), 1.0) ) print_error(" wrong CPXchgcoef [y11]");
            }
        }
    }
    
    // sum {i,j != i}^n [y(i,j,t)] = 1, for each t
    for (int t = 0; t < inst -> nnodes; t++) {
        int lastrow = CPXgetnumrows(env,lp);
        double rhs = 1.0;
        char sense = 'E';
        sprintf(cname[0], "constr2(%d)", t+1);
        if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) print_error(" wrong CPXnewrows [y1]");
        for (int i = 0; i < inst -> nnodes; i++) {
            for (int j = 0; j < inst -> nnodes; j++) {
                if (j == i) continue;
                if ( CPXchgcoef(env, lp, lastrow, y3pos(i, j, t, inst), 1.0) ) print_error(" wrong CPXchgcoef [y11]");
            }
        }
    }
    
    // sum {j,t >=2}^n [t * y(i,j,t)] - sum {k,t}^n [t * y(k,i,t)] = 1
    for (int i = 1; i < inst -> nnodes; i++) {
        int lastrow = CPXgetnumrows(env,lp);
        double rhs = 1.0;
        char sense = 'E';
        sprintf(cname[0], "constr3(%d)", i+1);
        if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) print_error(" wrong CPXnewrows [y1]");
        for ( int j = 0; j < inst -> nnodes; j++ ) {
            for ( int k = 1; k < inst -> nnodes; k++ ) {
                if ( CPXchgcoef(env, lp, lastrow, y3pos(i, j, k, inst), k) ) print_error(" wrong CPXchgcoef [y11]");
                if ( CPXchgcoef(env, lp, lastrow, y3pos(j, i, k, inst), - k) ) print_error(" wrong CPXchgcoef [y11]");
            }
        }
        for ( int j = 1; j < inst -> nnodes; j++ ) {
            if ( CPXchgcoef(env, lp, lastrow, y3pos(j, i, 0, inst), - 1) ) print_error(" wrong CPXchgcoef [y11]");
        }
    }
    
    // ********************************************************************** //
    // ************************** LAZY CONSTRAINTS ************************** //
    // ********************************************************************** //
    
    int izero = 0;
    int index[2];
    double value[2];
    double rhs = 0.0;
    char sense = 'E';
    
    // sum {i=1}^n [y(i,1,k)] = 0, for each k - {1}
    for (int i = 1; i < inst -> nnodes; i++) {
        for (int j = 0; j < inst -> nnodes; j++) {
            sprintf(cname[0], "constr1(%d,%d)", i+1, j+1);
            index[0] = y3pos(i, j, 0, inst);
            value[0] = 1.0;
            index[1] = y3pos(j, i, inst -> nnodes - 1, inst);
            value[1] = - 1.0;
            if ( CPXaddlazyconstraints(env, lp, 1, 2, &rhs, &sense, &izero, index, value, cname) ) print_error("wrong CPXlazyconstraints on 2-node SECs");
        }
    }
    
    // without lazy constraint
    
    /*
     // sum {i=1}^n [y(i,1,k)] = 0, for each k - {1}
     for (int i = 1; i < inst -> nnodes; i++) {
     for (int j = 0; j < inst -> nnodes; j++) {
     int lastrow = CPXgetnumrows(env,lp);
     double rhs = 0.0;
     char sense = 'E';
     sprintf(cname[0], "constr3(%d,%d)", i+1, j+1);
     if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) print_error(" wrong CPXnewrows [y1]");
     if ( CPXchgcoef(env, lp, lastrow, y3pos(i, j, 0, inst), 1.0) ) print_error(" wrong CPXchgcoef [y11]");
     if ( CPXchgcoef(env, lp, lastrow, y3pos(j, i, inst -> nnodes - 1, inst), - 1.0) ) print_error(" wrong CPXchgcoef [y11]");
     }
     } */
    
    // save the model
    if ( VERBOSE >= -100 ) CPXwriteprob(env, lp, "model.lp", NULL);
    
    free(cname[0]);
    free(cname);
    
}

/**
 Model 7: compact model T3
 
 @param inst instance of the struct "instance" for TSP problem.
 @param env CPLEX environment.
 @param lp CPLEX LP.
 */
void build_model_7(instance *inst, CPXENVptr env, CPXLPptr lp) {
    
    printf("Resolve instance \"%s\" with Compact Model: T3\n\n", inst -> input_file);
    
    // type: Binary
    char binary = 'B';
    // type: Continuous
    char continuous = 'C';
    
    char **cname = (char **) calloc(1, sizeof(char*));
    cname[0] = (char *) calloc(100, sizeof(char));
    
    // ********************************************************************** //
    // ******************************* VARIABLES **************************** //
    // ********************************************************************** //
    
    // add binary var x(i,j) such that 0 <= x(i,j) <= 1
    for ( int i = 0; i < inst->nnodes; i++ )
    {
        for ( int j = 0; j < inst->nnodes; j++ )
        {
            sprintf(cname[0], "x(%d,%d)", i+1,j+1);
            double obj = dist(i,j,inst); // cost == distance
            double lb = 0.0;
            double ub = ( i == j ) ? 0.0 : 1.0;
            if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname) ) print_error(" wrong CPXnewcols on x var.s");
            if ( CPXgetnumcols(env,lp)-1 != xpos_compact(i, j, inst) ) print_error(" wrong position for x var.s");
        }
    }
    
    // get the last column
    inst -> ystart = CPXgetnumcols(env, lp);
    
    // add continuous var y(i,j) such that 0 <= y(i,j) <= n
    for (int t = 0; t < inst -> nnodes; t++) {
        for ( int i = 0; i < inst -> nnodes; i++ ) {
            for ( int j = 0; j < inst -> nnodes; j++ ) {
                sprintf(cname[0], "y(%d,%d,%d)",t+1, i+1, j+1);
                double obj = 0.0;
                double lb = 0.0;
                double ub = inst -> nnodes - 1;
                if (j == 0 && t != (inst -> nnodes - 1)) ub = 0;
                if (i == 0 && t != 0) ub = 0;
                if (t == 0 && i != 0) ub = 0;
                if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &continuous, cname) ) print_error(" wrong CPXnewcols on y var.s");
                if ( CPXgetnumcols(env,lp) - 1 != y3pos(i, j, t, inst) ) print_error(" wrong position for y var.s");
            }
        }
    }
    
    // ********************************************************************** //
    // ***************************** CONSTRAINTS **************************** //
    // ********************************************************************** //
    
    int *cindex = (int *) calloc(inst -> nnodes, sizeof(int));
    double *cvalue = (double *) calloc(inst -> nnodes, sizeof(double));
    int begin[1];
    begin[0] = 0;
    
    // x-constraints (out- and in-degree at nodes)
    for ( int h = 0; h < inst -> nnodes; h++ )  // out-degree
    {
        double rhs = 1.0;
        char sense = 'E';
        sprintf(cname[0], "outdeg(%d)", h+1);
        for ( int i = 0; i < inst->nnodes; i++ )
        {
            cindex[i] = xpos_compact(i, h, inst);
            cvalue[i] = 1.0;
        }
        if (CPXaddrows(env, lp, 0, 1, inst -> nnodes, &rhs, &sense, begin, cindex, cvalue, NULL, cname)) print_error("Error in add row out-deg");
    }
    
    for ( int h = 0; h < inst->nnodes; h++ )   // in-degree
    {
        double rhs = 1.0;
        char sense = 'E';
        sprintf(cname[0], "indeg(%d)", h+1);
        for ( int i = 0; i < inst->nnodes; i++ )
        {
            cindex[i] = xpos_compact(h, i, inst);
            cvalue[i] = 1.0;
        }
        if (CPXaddrows(env, lp, 0, 1, inst -> nnodes, &rhs, &sense, begin, cindex, cvalue, NULL, cname)) print_error("Error in add row in-deg");
    }
    
    free(cindex);
    free(cvalue);
    
    // x(i,j) - sum {t}^n [y(i,j,t)] = 0
    for (int i = 0; i < inst -> nnodes; i++) {
        for (int j = 0; j < inst -> nnodes; j++) {
            if( i == j ) continue;
            int lastrow = CPXgetnumrows(env,lp);
            double rhs = 0.0;
            char sense = 'E';
            sprintf(cname[0], "constr(%d,%d)", i+1, j+1);
            if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) print_error(" wrong CPXnewrows [y1]");
            if ( CPXchgcoef(env, lp, lastrow, xpos_compact(i, j, inst), 1.0) ) print_error(" wrong CPXchgcoef [y11]");
            for ( int t = 0; t < inst -> nnodes; t++ ) {
                if ( CPXchgcoef(env, lp, lastrow, y3pos(i, j, t, inst), - 1.0) ) print_error(" wrong CPXchgcoef [y11]");
            }
        }
    }
    
    // sum {j}^n [y(1,j,1)] = 1
    int lastrow = CPXgetnumrows(env,lp);
    double rhs = 1.0;
    char sense = 'E';
    sprintf(cname[0], "constr(%d)", 1);
    if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) print_error(" wrong CPXnewrows [y1]");
    for ( int j = 0; j < inst -> nnodes; j++ ) {
        if ( CPXchgcoef(env, lp, lastrow, y3pos(0, j, 0, inst), 1.0) ) print_error(" wrong CPXchgcoef [y11]");
    }
    
    // sum {i}^n [y(i,1,n)] = 1
    lastrow = CPXgetnumrows(env,lp);
    sprintf(cname[0], "constr1(%d)", 1);
    if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) print_error(" wrong CPXnewrows [y1]");
    for ( int i = 0; i < inst -> nnodes; i++ ) {
        if ( CPXchgcoef(env, lp, lastrow, y3pos(i, 0, inst -> nnodes - 1, inst), 1.0) ) print_error(" wrong CPXchgcoef [y11]");
    }
    
    // sum {j}^n [y(i,j,t)] - sum {k}^n [y(k,i,t - 1)] = 0, for each i,t - {1}
    for ( int i = 1; i < inst -> nnodes; i++ ) {
        for ( int t = 1; t < inst -> nnodes; t++ ) {
            int lastrow = CPXgetnumrows(env,lp);
            double rhs = 0.0;
            char sense = 'E';
            sprintf(cname[0], "constr1(%d,%d)", i+1, t+1);
            if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) print_error(" wrong CPXnewrows [y1]");
            for ( int j = 0; j < inst -> nnodes; j++ ) {
                if ( CPXchgcoef(env, lp, lastrow, y3pos(i, j, t, inst), 1.0) ) print_error(" wrong CPXchgcoef [y11]");
            }
            for ( int k = 0; k < inst -> nnodes; k++ ) {
                if ( CPXchgcoef(env, lp, lastrow, y3pos(k, i, t - 1, inst), - 1.0) ) print_error(" wrong CPXchgcoef [y11]");
            }
        }
    }
    
    // ********************************************************************** //
    // ************************** LAZY CONSTRAINTS ************************** //
    // ********************************************************************** //
    
    int izero = 0;
    int index[2];
    double value[2];
    rhs = 0.0;
    sense = 'E';
    
    // sum {i=1}^n [y(i,1,k)] = 0, for each k - {1}
    for (int i = 1; i < inst -> nnodes; i++) {
        for (int j = 0; j < inst -> nnodes; j++) {
            sprintf(cname[0], "constr2(%d,%d)", i+1, j+1);
            index[0] = y3pos(i, j, 0, inst);
            value[0] = 1.0;
            index[1] = y3pos(j, i, inst -> nnodes - 1, inst);
            value[1] = - 1.0;
            if ( CPXaddlazyconstraints(env, lp, 1, 2, &rhs, &sense, &izero, index, value, cname) ) print_error("wrong CPXlazyconstraints on 2-node SECs");
        }
    }
    
    // save the model
    if ( VERBOSE >= -100 ) CPXwriteprob(env, lp, "model.lp", NULL);
    
    free(cname[0]);
    free(cname);
    
}


// ********************************************************************** //
// ************************* BUILD MODEL SOLUTION *********************** //
// ********************************************************************** //

/**
 Build succ() and comp() wrt xstar().

 @param xstar solution in CPLEX format.
 @param inst instance of the struct "instance" for TSP problem.
 @param succ solution as successors.
 @param comp number of the component associated to each node in the solution.
 @param ncomp number of components in the solution.
 */
void build_compact_sol(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp) {
    
    // only for debug
#if debug
    int *degree = (int *) calloc(inst->nnodes, sizeof(int));
    for ( int i = 0; i < inst->nnodes; i++ )
    {
        for ( int j = 0; j < inst->nnodes; j++ )
        {
            int k = xpos_compact(i,j,inst);
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
                if((i != j) && xstar[xpos_compact(i,j,inst)] > 0.5 && comp[j] == -1) {
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
