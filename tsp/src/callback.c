//
//  callback.c
//  cplex
//
//  Created by Matteo Moratello on 25/04/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//

#include "callback.h"

void print_error(const char *err);
// position
int xpos(int i, int j, instance *inst);
void build_sol(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp);

// ****************************** MODELS DEFINITION ****************************** //
//
// Lazy callback
// Add SEC constraints calling the CPEX callback in integer solutions
// User cut callback
// Add SEC constraints calling the CPEX callback in relaxation
// both are used in order to cut in efficient way the branching tree
//
// Generic callback
// It substitute the two previous callbacks with a single one (the callback does not deactivate the dinamic search algorithm)
//
// In both versions in the relaxation callbacks the Concorde's algorithms are used
// First check is the graph is connected, is affermative obtain all the cuts usign an efficient Flow algorithm.
// All found cuts are added to the model with the callback.
//
// ******************************************************************************* //


// ********************************* LAZY CALLBACK ******************************* //

// lazy callback: integer LP
static int CPXPUBLIC mylazycallback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p) {
    
    *useraction_p = CPX_CALLBACK_DEFAULT;
    instance* inst = (instance *) cbhandle;             // casting of cbhandle
    
    // get solution xstar
    double *xstar = (double*) malloc(inst -> ncols * sizeof(double));
    if ( CPXgetcallbacknodex(env, cbdata, wherefrom, xstar, 0, inst-> ncols - 1) ) return 1; // xstar = current x from CPLEX -- xstar starts from position 0
    
    // get some random information at the node (as an example)
    double objval = CPX_INFBOUND; CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &objval);
    int mythread = -1; CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_MY_THREAD_NUM, &mythread);
    double zbest; CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_BEST_INTEGER, &zbest);
    
    printf("Objval: %lf \nThread: %d \n", objval, mythread);
    
    //apply cut separator and possibly add violated cuts
    int ncuts = myseparation(inst, xstar, env, cbdata, wherefrom);
    printf("ncuts: %d\n\n", ncuts);
    
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
        
    }
    
    free(comp);
    free(succ);
    
    return (ncomp == 1? 0 : ncomp);
}

// User callback
static int CPXPUBLIC UserCutCallback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p) {
    
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
    printf("cutval: %lf \ncutcount: %d \n", cutval, cutcount);
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
            if (CPXcallbackrejectcandidate(context, 1, nnz, &rhs, &sense, start_indexes, index, value)) print_error("USER_separation: CPXcutcallbackadd error");
        }
        
        free(start_indexes);
        free(value);
        free(index);
        
        free(cname[0]);
        free(cname);
        
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






