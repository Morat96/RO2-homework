//
//  localbranching.c
//  cplex
//
//  Created by Matteo Moratello on 27/04/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//

#include "localbranching.h"

// ********************************************************************** //
// *************************** LOCAL BRANCHING ************************** //
// ********************************************************************** //


// ************************** MODEL DEFINITION ************************** //
//
// obj: compute an optimal solution within total time limit -> heuristic method
//
// Compute the solution within the internal timelimit (= total timelimit / 5)
// Update remaining time
// If the current solution doesn't improve too much (+10%) for two consecutive times
// change the rhs
// Repeat the sequence until the end of time
//
// Parameters to tune: timelimit and percentages
//
// ********************************************************************** //
void localbranching(instance *inst, CPXENVptr env, CPXLPptr lp) {
    
    printf("Resolve instance \"%s\" with Local Branching\n\n", inst -> input_file);
    
    double percentage = inst -> nnodes * 0.95;
    double remaining_time = inst -> timelimit;
    double time_x_cycle = (inst -> timelimit) / 5;
    double current_obj_val = CPX_INFBOUND;
    int cnt = 0;
    printf("Time x cycle: %lf seconds\n\n", time_x_cycle);
    double start = second();
    int cnt_cycle = 0;
    
    while (remaining_time >= 0) {
        
        // set internal timelimit
        if (remaining_time > time_x_cycle) {
            if (CPXsetdblparam(env, CPX_PARAM_TILIM, time_x_cycle)) print_error("error in change timelimit");
        }
        else {
            if (CPXsetdblparam(env, CPX_PARAM_TILIM, remaining_time)) print_error("error in change timelimit");
        }
        
        // take start time for computing total time used by CPXmipopt
        double t1 = second();
        
        // ** Find the best solution with CALLBACK METHOD within this time ** //
        CPXcallbacksetfunc(env, lp, CPX_CALLBACKCONTEXT_CANDIDATE | CPX_CALLBACKCONTEXT_RELAXATION, my_generic_callback, inst);
        int ncores = 1; CPXgetnumcores(env, &ncores);
        CPXsetintparam(env, CPX_PARAM_THREADS, ncores);
        
        inst -> ncols = CPXgetnumcols(env,lp);
        
        // compute the solution
        if (CPXmipopt(env,lp)) print_error("Error in find a solution");
        
        // take final time taken by CPXmipopt
        double t2 = second();
        // time up to now
        printf("Time: %lf\n", t2 - start);
        // time for computing the solution
        printf("Time for computing solution: %lf\n", t2 - t1);
        // remaining time
        remaining_time -= (t2-t1);
        printf("Remaining time : %lf\n", remaining_time);
        
        if( remaining_time <= 0) break;
        
        // Obtain the solution
        double objval;
        // value of objective function
        if (CPXgetobjval (env, lp, &objval)) print_error("Error in CPXgetobjval");
        printf("Objective function value: %.0f\n" , objval);
        
        // if the current solution doesn't improve too much (+10%) for two consecutive times
        // change the percentage of fixed variables
        if (objval > (current_obj_val - current_obj_val * 0.1)) cnt ++;
        printf("CNT: %d\n", cnt);
        
        if( cnt_cycle ) {
            if (percentage >= 0.1) percentage -= (percentage * 0.05);
            else percentage = 0.0;
            cnt = 0;
        }
        printf("Current rhs: %.0f\n\n", percentage);
        current_obj_val = objval;
        
        // number of variables of the problem
        int ncols = CPXgetnumcols(env, lp);
        
        // final value of variables
        double *xstar = (double *) calloc(ncols, sizeof(double));
        if (CPXgetx(env, lp, xstar, 0, ncols-1)) print_error("Error in CPXgetx");
        
        double *value_sol = (double *) calloc(inst -> nnodes, sizeof(double));
        int *index_sol = (int *) calloc(inst -> nnodes, sizeof(int));
        int ind_sol = 0;
        
        // reset the lower bound of variables to 0.0 (in the first take of first step is useless)
        int *ind = (int *) calloc(inst -> ncols, sizeof(int));
        char *lu = (char *) calloc(inst -> ncols, sizeof(char));
        double *db = (double *) calloc(inst -> ncols, sizeof(double));
        
        int rows = CPXgetnumrows(env, lp);
        
        // remove last row, which is the local branching constraint
        if (cnt_cycle++) CPXdelrows(env, lp, rows - 1, rows - 1);
        
        // save indices of variables of the solution
        for (int k = 0; k < ncols; k++) {
            if(xstar[k] > 0.5) {
                index_sol[ind_sol] = k;
                value_sol[ind_sol++] = 1.0;
            }
        }
        
        // ADD the new local branching constraint
        char **cname = (char **) calloc(1, sizeof(char*));
        cname[0] = (char *) calloc(100, sizeof(char));
        sprintf(cname[0], "local_branch");
        
        char sense = 'G';
        int begin[1];
        begin[0] = 0;
        
        if (CPXaddrows(env, lp, 0, 1, inst -> nnodes, &percentage, &sense, begin, index_sol, value_sol, NULL, cname)) print_error("Error in add rows");
        
        free(cname);
        free(value_sol);
        free(index_sol);
        free(ind);
        free(lu);
        free(db);
        
        free(xstar);
    } // repeat
}


// ************************** MODEL DEFINITION ************************** //
//
// obj: compute an optimal solution within total time limit -> heuristic method
// try to solve the problem dividing the total time in 5 time series
// --> time x step = total time / 5
// in each step add a constraint to the model --> sum {x(i,j) = 1} x(i,j) >= k * n
// if the obj_value[n] == obj_value[n - 1] --> go to the next step
// Each step has (possible) different value of k
// if the step ends earlier, go to the next step
//
// Parameters to tune: timelimit and k
//
// ********************************************************************** //
void localbranching_2nd_vers(instance *inst, CPXENVptr env, CPXLPptr lp) {
    
    // local time limit
    double time_x_cycle = (inst -> timelimit) / 5;
    printf("Time x cycle: %lf seconds\n\n", time_x_cycle);
    
    // change k of each step
    double percentage[5];
    percentage[0] = inst -> nnodes * 0.8;
    percentage[1] = inst -> nnodes * 0.8;
    percentage[2] = inst -> nnodes * 0.8;
    percentage[3] = inst -> nnodes * 0.9;
    percentage[4] = inst -> nnodes * 0.9;
    
    double rhs = inst -> nnodes * 0.9;
    
    // save the last two solution in order to compare the gap between them
    double solutions[2];
    double last_sol = 0.0;
    
    // 5 steps
    for (int s = 0; s < 5; s++) {
        
        // set local timelimit
        if (CPXsetdblparam(env, CPX_PARAM_TILIM, time_x_cycle)) print_error("error in change timelimit");
        
        // take start time of the step
        double t1 = second();
        
        printf("Iteration %d\n", s + 1);
        
        for (int t = 0; t < 3; t++) {
            
            // ** Find the best solution with CALLBACK METHOD within this time ** //
            
            // Second Version: generic callback
            CPXcallbacksetfunc(env, lp, CPX_CALLBACKCONTEXT_CANDIDATE | CPX_CALLBACKCONTEXT_RELAXATION, my_generic_callback, inst);
            
            int ncores = 1; CPXgetnumcores(env, &ncores);
            CPXsetintparam(env, CPX_PARAM_THREADS, ncores);
            
            inst -> ncols = CPXgetnumcols(env,lp);
            
            if (CPXmipopt(env,lp)) print_error("Error in find a solution");
            
            // take time taken by CPXmipopt
            double t2 = second();
            
            printf("Time: %lf/%lf\n", t2 - t1, time_x_cycle);
            if (((t2 - t1) > time_x_cycle) && s == 4) break;
            
            // compute the remaining time for the current step
            double remaining_time = time_x_cycle - (t2 - t1);
            if (remaining_time > 0) {
                // update the timelimit
                if (CPXsetdblparam(env, CPX_PARAM_TILIM, remaining_time)) print_error("error in change timelimit");
                printf("Remaining Time: %lf\n", remaining_time);
            }
            else {
                printf("Time limit reached! \n");
            }
            
            // Obtain the solution
            double objval;
            // value of objective function
            if (CPXgetobjval (env, lp, &objval)) print_error("Error in CPXgetobjval");
            printf("Objective function value: %lf\n" , objval);
            
            // if final step and solution not improve, exit
            if (last_sol == solutions[1] && s == 4) break;
            
            // obtain informations about last solution
            if (!( s == 0 && t == 0)) {
                last_sol = solutions[0];
                solutions[1] = objval;
            }
            solutions[0] = objval;
            
            // number of variables of the problem
            int ncols = CPXgetnumcols(env, lp);
            
            // final value of variables
            double *xstar = (double *) calloc(ncols, sizeof(double));
            if (CPXgetx(env, lp, xstar, 0, ncols-1)) print_error("Error in CPXgetx");
            
            int *index_sol = (int *) calloc(inst -> nnodes, sizeof(int));
            double *value_sol = (double *) calloc(inst -> nnodes, sizeof(double));
            int ind_sol = 0;
            
            int rows = CPXgetnumrows(env, lp);
            
            // remove last row, which is the local branching constraint
            if (!( s == 0 && t == 0) && !(s == 4 && t == 3)) {
                CPXdelrows(env, lp, rows - 1, rows - 1);
            }
            
            // save indices of variables of the solution
            for (int k = 0; k < ncols; k++) {
                if(xstar[k] > 0.5) {
                    index_sol[ind_sol] = k;
                    value_sol[ind_sol++] = 1.0;
                }
            }
            
            // set for each step the correct k
            if (s == 0) rhs = percentage[0];
            if (s == 1) rhs = percentage[1];
            if (s == 2) rhs = percentage[2];
            if (s == 3) rhs = percentage[3];
            if (s == 4) rhs = percentage[4];
            printf("Current rhs: %f\n\n", rhs);
            
            // ADD the new local branching constraint
            char **cname = (char **) calloc(1, sizeof(char*));
            cname[0] = (char *) calloc(100, sizeof(char));
            sprintf(cname[0], "local_branch");
            
            char sense = 'G';
            int begin[1];
            begin[0] = 0;
            
            if (!(s == 4 && t == 3)) {
                if (CPXaddrows(env, lp, 0, 1, inst -> nnodes, &rhs, &sense, begin, index_sol, value_sol, NULL, cname)) print_error("Error in add rows");
            }
            
            // if there is no improvment in the solution, go to the next step
            if (last_sol == solutions[1]) break;
            // if timelimit is reached, go to the next step
            if ((t2 - t1) > time_x_cycle) break;
            
            free(cname);
            free(value_sol);
            free(index_sol);
            free(xstar);
        }
    }
}
