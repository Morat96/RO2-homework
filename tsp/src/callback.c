//
//  callback.c
//  cplex
//
//  Created by Matteo Moratello on 25/04/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//

#include "callback.h"

// ****************************** MODELS DEFINITION ****************************** //
//
// Lazy callback
// Add SEC constraints calling the CPEX callback in integer solutions
// User cut callback
// Add SEC constraints calling the CPEX callback in relaxation
// both are used in order to cut in efficient way the branching tree
// Heuristic callback
// build and provide to CPLEX a TSP solution from an integer one
//
// Generic callback
// It substitute the three previous callbacks with a single one (the callback does not deactivate the dinamic search algorithm)
//
// In both versions in the relaxation callbacks the Concorde's algorithms are used
// First check is the graph is connected, is affermative obtain all the cuts usign an efficient Flow algorithm.
// All found cuts are added to the model with the callback.
//
// ******************************************************************************* //


// ********************************* LAZY CALLBACK ******************************* //

// lazy callback: integer LP
int CPXPUBLIC mylazycallback(CPXCENVptr env,
                             void *cbdata,
                             int wherefrom,
                             void *cbhandle,
                             int *useraction_p) {
    
    *useraction_p = CPX_CALLBACK_DEFAULT;
    instance* inst = (instance *) cbhandle;             // casting of cbhandle
    
    // get solution xstar
    double *xstar = (double*) malloc(inst -> ncols * sizeof(double));
    if ( CPXgetcallbacknodex(env, cbdata, wherefrom, xstar, 0, inst-> ncols - 1) ) return 1; // xstar = current x from CPLEX -- xstar starts from position 0
    
    // get some random information at the node (as an example)
    double objval = CPX_INFBOUND; CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &objval);
    int mythread = -1; CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_MY_THREAD_NUM, &mythread);
    double zbest; CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_BEST_INTEGER, &zbest);
    
    //printf("Objval: %lf \nThread: %d \n", objval, mythread);
    
    //apply cut separator and possibly add violated cuts
    int ncuts = myseparation(inst, xstar, env, cbdata, wherefrom);
    //printf("ncuts: %d\n\n", ncuts);
    
    free(xstar);
    
    if (ncuts > 1) *useraction_p = CPX_CALLBACK_SET;        // tell CPLEX that cuts have been created*/
    return 0;                                               // return 1 would mean error --> abort Cplex's execution
}

// add SEC constraints
int myseparation(instance *inst, double *xstar, CPXCENVptr env, void *cbdata, int wherefrom) {
    
    int *succ = (int *) calloc(inst->nnodes, sizeof(int));
    int *comp = (int *) calloc(inst->nnodes, sizeof(int));
    int ncomp = 9999;
    
    build_sol(xstar, inst, succ, comp, &ncomp);
    
    if (ncomp > 1) {
        
        char **cname = (char **) calloc(1, sizeof(char*));
        cname[0] = (char *) calloc(100, sizeof(char));
        
        // max (n * (n - 1) / 2) new constraints
        int *index = (int *) calloc((inst->nnodes * inst->nnodes - 1)/2, sizeof(int));
        double *value = (double *) calloc((inst->nnodes * inst->nnodes - 1)/2, sizeof(double));
        
        int nnz = 0;
        double rhs = 0.0;
        char sense = 'L';
        
        // for each component
        for(int i = 0; i < ncomp; i++) {
            // SEC of "n" iteration and "i+1" component
            sprintf(cname[0], "SEC(%d)", i+1);
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
            if ( CPXcutcallbackadd(env, cbdata, wherefrom, nnz, rhs, sense, index, value, 0) ) print_error("USER_separation: CPXcutcallbackadd error");
        }
        
        free(index);
        free(value);
        
        free(cname[0]);
        free(cname);
        
        // total number of nodes solved
        CPXINT node_count = 0; if (CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_NODE_COUNT, &node_count)) print_error("USER_separation: get info node count error");
        
        // provide to MIP a TSP integer solution only in nodes within depth nine of the branching tree
        // or every 10 nodes
        if (node_count < 512 || (node_count % 10 == 0)) {
            
            // save the complete graph from components found by CPLEX
            // each thread has a specific solution
            // flag[i] = 1 -> thread i has a solution available
            int mythread = -1; if (CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_MY_THREAD_NUM, &mythread)) print_error("USER_separation: get info thread error");
            complete_cycle(inst, succ, comp, &ncomp);
            for (int i = 0; i < inst -> nnodes; i++) inst -> sol_thread[mythread][i] = succ[i];
            inst -> flag[mythread] = 1;
        }
    }
    
    free(comp);
    free(succ);
    
    return (ncomp == 1? 0 : ncomp);
}

// User callback
int CPXPUBLIC UserCutCallback(CPXCENVptr env,
                              void *cbdata,
                              int wherefrom,
                              void *cbhandle,
                              int *useraction_p) {
    
    *useraction_p = CPX_CALLBACK_DEFAULT;
    instance* inst = (instance *) cbhandle;                                                  // casting of cbhandle
    
    // get solution xstar
    double *xstar = (double*) malloc(inst -> ncols * sizeof(double));
    if ( CPXgetcallbacknodex(env, cbdata, wherefrom, xstar, 0, inst-> ncols - 1) ) return 1; // xstar = current x from CPLEX -- xstar starts from position 0
    
    int nedge = inst -> nnodes * ( inst -> nnodes - 1) / 2;
    int *elist = (int*) malloc ((nedge * 2) * sizeof(int));
    int loader = 0;
    for (int i = 0; i < inst -> nnodes; ++i) {
        for (int j = i + 1; j < inst -> nnodes; ++j) {
            elist[loader++] = i;
            elist[loader++] = j;
        }
    }
    int ncomp = 0;
    int *comps = (int*) malloc (inst -> nnodes * sizeof(int));
    int *compscount = (int*) malloc (inst -> nnodes * sizeof(int));
    if ( CCcut_connect_components(inst -> nnodes, nedge, elist, xstar, &ncomp, &compscount, &comps)) print_error( "Error during concorde connect comps algortihm!") ;
    
    //printf("components found: %d\n", ncomp);
    const float EPSILON = 0.1;
    
    input in;
    in.inst = inst;
    in.env = env;
    in.cbdata = cbdata;
    in.wherefrom = wherefrom;
    in.useraction_p = useraction_p;
    
    if (ncomp == 1) {
        if ( CCcut_violated_cuts (inst -> nnodes, nedge, elist, xstar, 2.0 - EPSILON, doit_fn_concorde, (void*)&in)) print_error("error in CCcut_violated_cuts");
        *useraction_p = CPX_CALLBACK_SET;
    }
    
    free(elist);
    free(comps);
    free(compscount);
    free(xstar);
    
    return 0;                                                                           // return 1 would mean error --> abort Cplex's execution
}

// add SEC in continuous relaxation solutions
int doit_fn_concorde(double cutval, int cutcount, int *cut , void *inParam) {
    //printf("cutval: %lf \ncutcount: %d \n", cutval, cutcount);
    input* inputVal = (input *) inParam;
    int *index = (int *) calloc(cutcount, sizeof(int));
    double *value = (double *) calloc(cutcount, sizeof(double));
    char sense = 'L';
    for (int i = 0; i < cutcount; i++) {
        index[i] = cut[i];
        value[i] = 1.0;
    }
    if ( CPXcutcallbackadd(inputVal -> env, inputVal -> cbdata, inputVal -> wherefrom, cutcount, cutcount - 1, sense, index, value, 0) ) print_error("USER_separation: CPXcutcallbackadd error");
    
    free(value);
    free(index);
    return 0;
}

// heuristic callback
// compute a heuristic solution from an integer solution of CPLEX
int CPXPUBLIC myheuristic (CPXCENVptr env,
                           void       *cbdata,
                           int        wherefrom,
                           void       *cbhandle,
                           double     *objval_p,
                           double     *x,
                           int        *checkfeas_p,
                           int        *useraction_p)
{
    
    int mythread = -1; CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_MY_THREAD_NUM, &mythread);
    instance* inst = (instance *) cbhandle;
    CPXINT node_count = 0; if (CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_NODE_COUNT, &node_count)) print_error("USER_separation: get info node count error");
    *useraction_p = CPX_CALLBACK_DEFAULT;
    
    if (inst -> flag[mythread]) {
        
        double objval = 0;
        
        for (int i = 0; i < inst -> ncols; i++) x[i] = 0.0;
        
        for (int i = 0; i < inst -> nnodes; i++) {
            int pos = xpos(i, inst -> sol_thread[mythread][i], inst);
            x[pos] = 1.0;
            objval += dist(i, inst -> sol_thread[mythread][i], inst);
        }
        
        printf("Obj value: %lf, Thread: %d, Node: %d\n", objval, mythread, node_count);
        inst -> flag[mythread] = 0;
        *objval_p = objval;
        *checkfeas_p = 1;
        *useraction_p = CPX_CALLBACK_SET;
    }
    
    return 0;
}

// ********************************************************************************* //

// ******************************** GENERIC CALLBACK ******************************* //

// generic callback implementation
int my_generic_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *user) {
    
    instance* inst = (instance *) user;
    
    // get solution xstar
    double *xstar = (double*) malloc(inst->ncols * sizeof(double));
    double objval = CPX_INFBOUND + 0.0;
    
    double mythread = -1; CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_THREADID, &mythread);
    double zbest; CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_BEST_SOL, &zbest);
    
    switch (contextid) {
        case CPX_CALLBACKCONTEXT_CANDIDATE: {
            
            if ( CPXcallbackgetcandidatepoint(context, xstar, 0, inst -> ncols - 1, &objval)) return 1;
            //printf("obj val (Candidate): %lf\n", objval);
            int ncuts;
            ncuts = my_separation(inst, xstar, context);
            //printf("ncuts: %d\n\n", ncuts);
            
            // compute a heuristic solution from an integer solution of CPLEX
            int mythread = -1; if (CPXcallbackgetinfoint(context, CPXCALLBACKINFO_THREADID, &mythread)) print_error("Error in generic callback: get info thread id error");
            CPXINT node_count = 0; if (CPXcallbackgetinfoint(context, CPXCALLBACKINFO_NODECOUNT, &node_count)) print_error("Error in generic callback: get info node count error");
            
            if (inst -> flag[mythread]) {
                
                double objval = 0;
                double* x = (double*) calloc(inst -> ncols, sizeof(double));
                int* ind = (int*) calloc(inst -> ncols, sizeof(int));
                for (int i = 0; i < inst -> ncols; i++) {
                    x[i] = 0.0;
                    ind[i] = i;
                }
                
                for (int i = 0; i < inst -> nnodes; i++) {
                    int pos = xpos(i, inst -> sol_thread[mythread][i], inst);
                    x[pos] = 1.0;
                    objval += dist(i, inst -> sol_thread[mythread][i], inst);
                }
                
                printf("Obj value: %lf, Thread: %d, Node: %d\n", objval, mythread, node_count);
                if (CPXcallbackpostheursoln(context, inst -> ncols, ind, x, objval, CPXCALLBACKSOLUTION_CHECKFEAS)) print_error("Error in generic callback: Heuristic");
                
                inst -> flag[mythread] = 0;
                free(x);
                free(ind);
            }
            break;
        }
        case CPX_CALLBACKCONTEXT_RELAXATION: {
            
            if (CPXcallbackgetrelaxationpoint(context, xstar, 0, inst -> ncols - 1, &objval)) return 1;
            
            //printf("obj val (Relaxation): %lf\n", objval);
            int nedge = inst -> nnodes * ( inst -> nnodes - 1) / 2;
            int *elist = (int*) malloc ((nedge * 2) * sizeof(int));
            int loader = 0;
            for (int i = 0; i < inst -> nnodes; ++i) {
                for (int j = i + 1; j < inst -> nnodes; ++j) {
                    elist[loader++] = i;
                    elist[loader++] = j;
                }
            }
            int ncomp = 0;
            int *comps = (int*) malloc (inst -> nnodes * sizeof(int));
            int *compscount = (int*) malloc (inst -> nnodes * sizeof(int));
            if ( CCcut_connect_components(inst -> nnodes, nedge, elist, xstar, &ncomp, &compscount, &comps)) print_error( "Error during concorde connect comps algortihm!") ;
            
            //printf("components found: %d\n", ncomp);
            const float EPSILON = 0.1;
            
            if (ncomp == 1) {
                if ( CCcut_violated_cuts (inst -> nnodes, nedge, elist, xstar, 2.0 - EPSILON, doit_fn_concorde_gen, (void*)context)) print_error("error in CCcut_violated_cuts");
            }
            
            free(elist);
            free(comps);
            free(compscount);
            
            break;
        }
    }
    
    free(xstar);
    
    return 0;
}

// add SEC constraints
int my_separation(instance *inst, double *xstar, CPXCALLBACKCONTEXTptr context) {
    
    int *succ = (int *) calloc(inst->nnodes, sizeof(int));
    int *comp = (int *) calloc(inst->nnodes, sizeof(int));
    int ncomp = 9999;
    
    build_sol(xstar, inst, succ, comp, &ncomp);
    
    if (ncomp > 1) {
        
        char **cname = (char **) calloc(1, sizeof(char*));
        cname[0] = (char *) calloc(100, sizeof(char));
        
        // max (n * (n - 1) / 2) new constraints
        int *index = (int *) calloc((inst->nnodes * inst->nnodes - 1)/2, sizeof(int));
        double *value = (double *) calloc((inst->nnodes * inst->nnodes - 1)/2, sizeof(double));
        int *start_indexes = (int *) calloc((inst->nnodes * inst->nnodes - 1)/2, sizeof(int));
        
        int nnz = 0;
        double rhs = 0.0;
        char sense = 'L';
        
        // for each component
        for(int i = 0; i < ncomp; i++) {
            // SEC of "n" iteration and "i+1" component
            sprintf(cname[0], "SEC(%d)", i+1);
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
            if (CPXcallbackrejectcandidate(context, 1, nnz, &rhs, &sense, start_indexes, index, value)) print_error("Generic callback - separation: CPXcutcallbackadd error");
        }
        
        free(start_indexes);
        free(value);
        free(index);
        
        free(cname[0]);
        free(cname);
        
        // total number of nodes solved
        CPXINT node_count = 0; if (CPXcallbackgetinfoint(context, CPXCALLBACKINFO_NODECOUNT, &node_count)) print_error("Generic callback - separation: get info node count error");
        
        // provide to MIP a TSP integer solution only in nodes within depth nine of the branching tree
        // or every 10 nodes
        if (node_count < 512 || (node_count % 10 == 0)) {
            
            // save the complete graph from components found by CPLEX
            // each thread has a specific solution
            // flag[i] = 1 -> thread i has a solution available
            int mythread = -1; if (CPXcallbackgetinfoint(context, CPXCALLBACKINFO_THREADID, &mythread)) print_error("Generic callback - separation: get info thread error");
            complete_cycle(inst, succ, comp, &ncomp);
            for (int i = 0; i < inst -> nnodes; i++) inst -> sol_thread[mythread][i] = succ[i];
            inst -> flag[mythread] = 1;
        }
    }
    
    free(succ);
    free(comp);
    
    return (ncomp == 1 ? 0 : ncomp);
}

// add SEC in continuous relaxation solutions
int doit_fn_concorde_gen(double cutval, int cutcount, int *cut , void *inParam) {
    
    //printf("Relaxation \ncutval: %lf \ncutcount: %d \n\n", cutval, cutcount);
    CPXCALLBACKCONTEXTptr context = (CPXCALLBACKCONTEXTptr) inParam;
    
    int *index = (int *) calloc(cutcount, sizeof(int));
    double *value = (double *) calloc(cutcount, sizeof(double));
    int *start_indexes = (int *) calloc(cutcount, sizeof(int));
    char sense = 'L';
    
    for (int i = 0; i < cutcount; i++) {
        index[i] = cut[i];
        value[i] = 1.0;
    }
    double rhs = cutcount - 1;
    int type_purge = CPX_USECUT_FILTER;
    CPXcallbackaddusercuts(context, 1, cutcount, &rhs, &sense, start_indexes, index, value, &type_purge, start_indexes);
    
    free(start_indexes);
    free(value);
    free(index);
    
    return 0;
}

// ********************************************************************************* //

// ****************************** MERGE COMP ALGORITHM ***************************** //

// build a cycle from the distinct components found by CPLEX
// merge components based on distances between nodes of different components
// until the solution has one component -> feasible solution
void complete_cycle(instance *inst, int *succ, int *comp, int *ncomp) {
    
    int first = 0;
    int second = 0;
    int succ_first = 0;
    int succ_second = 0;
    int comp_f = 0;
    int comp_s = 0;
    int* inv = (int *) calloc(inst -> nnodes, sizeof(int));
    int* index = (int *) calloc(inst -> nnodes, sizeof(int));
    
    while((*ncomp) != 1) {
        
        double min_function = INT_MAX;
        double curr_min_function = INT_MAX;
        int flag = 0;
        
        for (int i = 0; i < inst -> nnodes; i++) {
            for (int j = i + 1; j < inst -> nnodes; j++) {
                if (comp[i] != comp[j]) {
                    // Δ(a,b)
                    curr_min_function = dist(i, succ[j], inst) + dist(j, succ[i], inst) - dist(i, succ[i], inst) - dist(j, succ[j], inst);
                    
                    if (curr_min_function < min_function) {
                        comp_f = comp[i];
                        comp_s = comp[j];
                        first = i;                          // a
                        second = j;                         // b
                        succ_first = succ[i];               // a'
                        succ_second = succ[j];              // b'
                        min_function = curr_min_function;
                        flag = 0;
                    }
                    // Δ'(a,b)
                    curr_min_function = dist(i, j, inst) + dist(succ[j], succ[i], inst) - dist(i, succ[i], inst) - dist(j, succ[j], inst);
                    
                    if (curr_min_function < min_function) {
                        comp_f = comp[i];
                        comp_s = comp[j];
                        first = i;                          // a
                        second = j;                         // b
                        succ_first = succ[i];               // a'
                        succ_second = succ[j];              // b'
                        min_function = curr_min_function;
                        flag = 1;
                    }
                }
            }
        }
        
        if (flag) {
            
            // reverse the order of edges of the second component
            int cnt = 0;
            for (int i = 0; i < inst -> nnodes; i++)
                if (comp[i] == comp_s) cnt++;
            
            int cnt2 = 0;
            for (int i = 0; i < inst -> nnodes; i++) {
                if (comp[i] == comp_s) {
                    int val = succ[i];
                    for (int j = 0; j < cnt - 2; j++) val = succ[val];
                    
                    inv[cnt2] = val;
                    index[cnt2++] = i;
                }
            }
            
            for (int i = 0; i < cnt; i++) succ[index[i]] = inv[i];
            
            // update
            succ[first] = second;                               // (a,a') -> (a,b)
            succ[succ_second] = succ_first;                     // (b,b') -> (b',a')
            
        }
        else {
            // update
            succ[first] = succ_second;                          // (a,a') -> (a,b')
            succ[second] = succ_first;                          // (b,b') -> (b,a')
        }
        
        // update components
        for (int i = 0; i < inst -> nnodes; i++) if (comp[i] == comp_s) comp[i] = comp_f;
        (*ncomp) --;
    }
    
    free(inv);
    free(index);
}




