//
//  hardfixing.c
//  cplex
//
//  Created by Matteo Moratello on 26/04/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//

#include "hardfixing.h"

// ********************************************************************** //
// ***************************** HARD FIXING **************************** //
// ********************************************************************** //


// ************************** MODEL DEFINITION ************************** //
//
// obj: compute an optimal solution within total time limit -> heuristic method
//
// Compute the solution within the internal timelimit (= total timelimit / 5)
// Update remaining time
// If the current solution doesn't improve too much (+10%) for two consecutive times
// change the percentage of fixed variables
// Repeat the sequence until the end of time
//
// Parameters to tune: timelimit and percentages
//
// ********************************************************************** //
void hardfixing(instance *inst, CPXENVptr env, CPXLPptr lp) {
    
    printf("Resolve instance \"%s\" with Hard Fixing\n\n", inst -> input_file);
    
    float percentage = 0.9;
    double remaining_time = inst -> timelimit;
    double time_x_cycle = (inst -> timelimit) / 5;
    double current_obj_val = CPX_INFBOUND;
    int cnt = 0;
    printf("Time x cycle: %lf seconds\n\n", time_x_cycle);
    double start = second();
    
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
        printf("Objective function value: %lf\n" , objval);
        
        // if the current solution doesn't improve too much (+10%) for two consecutive times
        // change the percentage of fixed variables
        if (objval > (current_obj_val - current_obj_val * 0.1)) cnt ++;
        printf("CNT: %d\n", cnt);
        if (cnt == 2) {
            if (percentage >= 0.15) percentage -= 0.15;
            else percentage = 0.0;
            cnt = 0;
        }
        current_obj_val = objval;
        
        // number of variables of the problem
        int ncols = CPXgetnumcols(env, lp);
        
        // final value of variables
        double *xstar = (double *) calloc(ncols, sizeof(double));
        if (CPXgetx(env, lp, xstar, 0, ncols-1)) print_error("Error in CPXgetx");
        
        int *index_sol = (int *) calloc(inst -> nnodes, sizeof(int));
        int ind_sol = 0;
        
        // reset the lower bound of variables to 0.0 (in the first take of first step is useless)
        int *ind = (int *) calloc(inst -> ncols, sizeof(int));
        char *lu = (char *) calloc(inst -> ncols, sizeof(char));
        double *db = (double *) calloc(inst -> ncols, sizeof(double));
            
        for (int j = 0; j < inst -> ncols; j++) {
            ind[j] = j;
            lu[j] = 'L';
            db[j] = 0.0;
        }
            
        CPXchgbds(env, lp, inst -> ncols, ind, lu, db);
        
        // save indices of variables of the solution
        for (int k = 0; k < ncols; k++) {
            if(xstar[k] > 0.5) index_sol[ind_sol++] = k;
        }
        
        // randomize seed
        srand((unsigned int)time(NULL));
        
        printf("Current percentage of fixed variables: %.0f%%\n\n", percentage * 100);
        
        // compute the number of variables to be fixed
        float rv = 0.0;
        int count = 0;
        for (int k = 0; k < inst -> nnodes; k++) {
            rv = ((float)rand()/(float)(RAND_MAX));
            if (rv <= percentage) {
                ind[count] = index_sol[k];
                lu[count] = 'L';
                db[count++] = 1.0;
            }
        }
        
        CPXchgbds(env, lp, count, ind, lu, db);
        
        free(ind);
        free(lu);
        free(db);
        free(index_sol);
        
        free(xstar);
    } // repeat
}

// ************************** MODEL DEFINITION ************************** //
//
// obj: compute an optimal solution within total time limit -> heuristic method
// try to solve the problem dividing the total time in 5 time series
// --> time x step = total time / 5
// in each step fix a percentage of solution variables, according to the following scheme:
// —> step 1: 80 % -> three take (random pattern) max
// —> step 2: 80 % -> three take (random pattern) max
// —> step 3: 60 % -> three take (random pattern) max
// —> step 4: 40 % -> three take (random pattern) max
// —> step 5: 20 % -> three take (random pattern) max
// for each step, try with at most three different random pattern (take)
// if the step ends earlier, go to the next step
//
// Parameters to tune: timelimit and percentages
//
// ********************************************************************** //
void hardfixing_2nd_vers(instance *inst, CPXENVptr env, CPXLPptr lp) {
    
    // change percentages of each step
    float percentage[5];
    percentage[0] = 0.8;
    percentage[1] = 0.8;
    percentage[2] = 0.6;
    percentage[3] = 0.4;
    percentage[4] = 0.2;

    // steps
    int total_take_counter = 0;

    // take the time for each step
    double time_x_cycle = (inst -> timelimit) / 5;
    printf("Time x cycle: %lf seconds\n\n", time_x_cycle);
    float eps = 5.0;

    // 5 steps
    for (int s = 0; s < 5; s ++) {
        
        // initial time
        if (CPXsetdblparam(env, CPX_PARAM_TILIM, time_x_cycle)) print_error("error in change timelimit");
        
            // take start time of the step
            double t1 = second();
        
            printf("Iteration %d\n", s + 1);
        
            // three take for each step
            for (int t = 0; t < 3; t ++) {
                
                // ** Find the best solution with CALLBACK METHOD within this time ** //
                
                // Second Version: generic callback
                CPXcallbacksetfunc(env, lp, CPX_CALLBACKCONTEXT_CANDIDATE | CPX_CALLBACKCONTEXT_RELAXATION, my_generic_callback, inst);
                
                int ncores = 1; CPXgetnumcores(env, &ncores);
                CPXsetintparam(env, CPX_PARAM_THREADS, ncores);
                
                inst -> ncols = CPXgetnumcols(env,lp);
                
                if (CPXmipopt(env,lp)) print_error("Error in find a solution");
                
                // take time taken by CPXmipopt
                double t2 = second();
                
                total_take_counter ++;
                printf("Step: %d \n", total_take_counter);
                
                printf("Time: %lf/%lf\n", t2 - t1, time_x_cycle);
                if (((t2 - t1) > (time_x_cycle - eps)) && s == 4) break;
                
                // remaining time in this step
                double remaining_time = time_x_cycle - (t2 - t1);
                if (remaining_time > 0) {
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
                
                // number of variables of the problem
                int ncols = CPXgetnumcols(env, lp);
                
                // final value of variables
                double *xstar = (double *) calloc(ncols, sizeof(double));
                if (CPXgetx(env, lp, xstar, 0, ncols-1)) print_error("Error in CPXgetx");
                
                int *index_sol = (int *) calloc(inst -> nnodes, sizeof(int));
                int ind_sol = 0;
                
                // reset the lower bound of variables to 0.0 (in the first take of first step is useless)
                if (!(s == 4 && t == 2)) {
                    
                    int *ind = (int *) calloc(inst -> ncols, sizeof(int));
                    char *lu = (char *) calloc(inst -> ncols, sizeof(char));
                    double *db = (double *) calloc(inst -> ncols, sizeof(double));
                    
                    for (int j = 0; j < inst -> ncols; j++) {
                        ind[j] = j;
                        lu[j] = 'L';
                        db[j] = 0.0;
                    }
                    
                    CPXchgbds(env, lp, inst -> ncols, ind, lu, db);
                    
                    free(ind);
                    free(lu);
                    free(db);
                }
                
                // save indices of variables of the solution
                for (int k = 0; k < ncols; k++) {
                    if(xstar[k] > 0.5) index_sol[ind_sol++] = k;
                }
                
                // randomize seed
                srand((unsigned int)time(NULL));
                
                // set the percentage of variables to be fixed
                float perc = 0.8;
                if (s == 0) perc = percentage[0];
                if (s == 1) perc = percentage[1];
                if (s == 2) perc = percentage[2];
                if (s == 3) perc = percentage[3];
                if (s == 4) perc = percentage[4];
                
                // compute the number of variables to be fixed
                int cnt = (inst -> nnodes) * perc;
                printf("Current percentage of fixed variables: %f%%\n\n", perc * 100);
                
                int *ind = (int *) calloc(cnt, sizeof(int));
                char *lu = (char *) calloc(cnt, sizeof(char));
                double *db = (double *) calloc(cnt, sizeof(double));
                int r = 0;
                
                // fix variables applying lower bound equal to 1.0
                if (!(s == 4 && t == 2)) {
                    for (int j = 0; j < cnt; j++) {
                        r = rand() % (inst -> nnodes);
                        ind[j] = index_sol[r];
                        lu[j] = 'L';
                        db[j] = 1.0;
                    }
                    
                    CPXchgbds(env, lp, cnt, ind, lu, db);
                }
                
                free(ind);
                free(lu);
                free(db);
                
                // repeat
                free(xstar);
                
                // if time exceeded go to the next step
                if ((t2 - t1) > (time_x_cycle - eps)) break;
            }
    }
}
