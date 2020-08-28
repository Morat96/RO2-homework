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
// Add SECs calling the CPEX callback in integer solutions
// User cut callback
// Add SECs calling the CPEX callback in relaxation
// both are used in order to cut in efficient way the branching tree
// Heuristic callback
// build and provide to CPLEX a TSP solution from an integer one
//
// Generic callback
// It substitute the three previous callbacks with a single one (the callback does not deactivate the dinamic search algorithm)
//
// In both versions in the relaxation callbacks the Concorde's algorithms are used
// First check is the graph is connected, is affermative obtain all the cuts using an efficient Flow algorithm.
// All found cuts are added to the model with the callback.
//
// ******************************************************************************* //


// ********************************* LAZY CALLBACK ******************************* //

/**
 Lazy callback: produce a cut for the LP TSP problem in integer solutions.

 @param env Pointer to the CPLEX environment.
 @param cbdata Pointer to pass to functions that obtain callback-specific information.
 @param wherefrom An integer value reporting where in the optimization this function was called.
 @param cbhandle Pointer to private user data.
 @param useraction_p Pointer to an integer specifying the action for CPLEX to take at the completion of the user callback.
 @return The routine returns 0 (zero) if successful and nonzero if an error occurs.
 */
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

/**
 Add SEC constraints to the problem.

 @param inst instance of the struct "instance" for TSP problem.
 @param xstar TSP solution in CPLEX format.
 @param env Pointer to the CPLEX environment.
 @param cbdata Pointer to pass to functions that obtain callback-specific information.
 @param wherefrom An integer value reporting where in the optimization this function was called.
 @return The routine returns 0 (zero) if successful and nonzero if an error occurs.
 */
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
        //CPXINT node_count = 0; if (CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_NODE_COUNT, &node_count)) print_error("USER_separation: get info node count error");
        
        CPXINT node_depth = 0;
        if (CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_DEPTH, &node_depth)) print_error("USER_separation: get info node depth error");
        
        // provide to MIP a TSP integer solution only in nodes within depth 10 of the branching tree
        if (node_depth <= 10) {
            
            // save the complete graph from components found by CPLEX
            // each thread has a specific solution
            // flag[i] = 1 -> thread i has a solution available
            int mythread = -1; if (CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_MY_THREAD_NUM, &mythread)) print_error("USER_separation: get info thread error");
            // merge components of the solution with the complete tour algorithm
            complete_cycle(inst, succ, comp, &ncomp);
            for (int i = 0; i < inst -> nnodes; i++) inst -> sol_thread[mythread][i] = succ[i];
            inst -> flag[mythread] = 1;
        }
    }
    
    free(comp);
    free(succ);
    
    return (ncomp == 1? 0 : ncomp);
}

/**
 User callback: produce cuts for the TSP in fractional solutions.

 @param env Pointer to the CPLEX environment.
 @param cbdata Pointer to pass to functions that obtain callback-specific information.
 @param wherefrom An integer value reporting where in the optimization this function was called.
 @param cbhandle Pointer to private user data.
 @param useraction_p Pointer to an integer specifying the action for CPLEX to take at the completion of the user callback.
 @return The routine returns 0 (zero) if successful and nonzero if an error occurs.
 */
int CPXPUBLIC UserCutCallback(CPXCENVptr env,
                              void *cbdata,
                              int wherefrom,
                              void *cbhandle,
                              int *useraction_p) {
    
    *useraction_p = CPX_CALLBACK_DEFAULT;
    instance* inst = (instance *) cbhandle;                                                  // casting of cbhandle
    
    // total number of nodes solved
    CPXINT node_count = 0; if (CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_NODE_COUNT, &node_count)) print_error("USER_separation: get info node count error");
    
    CPXINT node_depth = 0;
    if (CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_DEPTH, &node_depth)) print_error("USER_separation: get info node depth error");
    //printf("Node depth: %d\n", node_depth);
    
    // add cuts in fractional solution only when the depth of the branching tree is less than or equal to 10
    if (node_depth <= 10) {
    
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
        
        //printf("Components found: %d\n", ncomp);
        const float EPSILON = 0.1;
        
        input in;
        in.inst = inst;
        in.env = env;
        in.cbdata = cbdata;
        in.wherefrom = wherefrom;
        in.useraction_p = useraction_p;
        
        if(ncomp > 1)
        {
            // Add constraints on connected components
            separationMultiple(inst, ncomp, compscount, comps, env, cbdata, wherefrom);
            *useraction_p = CPX_CALLBACK_SET;         // tell CPLEX that cuts have been created
        }
        if (ncomp == 1) {
            if ( CCcut_violated_cuts (inst -> nnodes, nedge, elist, xstar, 2.0 - EPSILON, doit_fn_concorde, (void*)&in)) print_error("error in CCcut_violated_cuts");
            *useraction_p = CPX_CALLBACK_SET;
        }
        
        free(elist);
        free(comps);
        free(compscount);
        free(xstar);
    }
    
    return 0;                                                                           // return 1 would mean error --> abort Cplex's execution
}

/**
 Add SEC in continuous relaxation solutions.
 
 @param cutval value of size of the cut.
 @param cutcount number of nodes in the cut.
 @param cut nodes belonging to the cut.
 @param inParam Pointer to private user data.
 @return The routine returns 0 (zero) if successful and nonzero if an error occurs.
 */
int doit_fn_concorde(double cutval, int cutcount, int *cut , void *inParam) {
    
    //printf("cutval: %lf \ncutcount: %d \n", cutval, cutcount);
    int i,j;
    input* inputVal = (input *) inParam;
    
    int num_x_var = (inputVal -> inst -> nnodes - 1) * inputVal -> inst -> nnodes / 2;     // number of x variables

    int nzcnt = 0;                    // number of non-zero variables in the constraint
    char sense = 'L';
    double rmatval[num_x_var];        // coefficients of the non-zero variables
    int rmatind[num_x_var];           // position of the variables to set (in terms of columns)
    int rhs = cutcount - 1;
    for(int i = 0; i < num_x_var; i++)
    {
        rmatval[i] = 0;
    }
    
    for(i = 0; i < cutcount; ++i)
    {
        for(j = 0; j < cutcount; ++j)
        {
            if(cut[i] >= cut[j])
            {
                continue;
            }
            else
            {
                rmatind[nzcnt] = xpos(cut[i], cut[j], inputVal->inst);
                rmatval[nzcnt] = 1.0;
                nzcnt++;
            }
        }
    }
    
    if(CPXcutcallbackadd(inputVal -> env, inputVal -> cbdata, inputVal -> wherefrom, nzcnt, rhs, sense, rmatind, rmatval, CPX_USECUT_FORCE)) {
        print_error("USER_separation: CPXcutcallbackadd error single component");
    }
    if(VERBOSE > 0) {
        printf("Cut added with one connected component\n");
    }

    return 0;
}

int separationMultiple(instance *inst, int ncomp, int *compscount, int *comps, CPXCENVptr env, void *cbdata, int wherefrom)
{
    int num_x_var = (inst -> nnodes - 1) * inst -> nnodes / 2;     // number of x variables
    int offset = 0;
    int added = 0;
    for(int i = 0; i < ncomp; i++) {
        
        int nzcnt = 0;                    // number of non-zero variables in the constraint
        char sense = 'L';
        double rmatval[num_x_var];        // coefficients of the non-zero variables
        int rmatind[num_x_var];           // position of the variables to set (in terms of columns)
        
        for(int i = 0; i < num_x_var; i++)
        {
            rmatval[i] = 0;
        }
        
        int comp_vertexes = compscount[i];
        double rhs = comp_vertexes - 1.0;
        
        int j,k;
        for(j = offset; j < offset + comp_vertexes; j++)
        {
            for(k = offset; k < offset + comp_vertexes; k++)
            {
                if((comps[j] == comps[k]) || (comps[j] > comps[k]))
                {
                    continue;
                }
                else
                {
                    rmatind[nzcnt] = xpos(comps[j], comps[k], inst);
                    rmatval[nzcnt] = 1.0;
                    nzcnt++;
                }
            }

        }

        offset += comp_vertexes;
        if(CPXcutcallbackadd(env, cbdata, wherefrom, nzcnt, rhs, sense, rmatind, rmatval, CPX_USECUT_FORCE)) {
            print_error("USER_separation: CPXcutcallbackadd error user callback");
        }
        added++;
    }
    if(VERBOSE > 0) {
        printf("Cut added with multiple connected component\n");
    }
    
    return ncomp;
}

/**
 Heuristic callback.
 Compute a heuristic solution from an integer solution of CPLEX.

 @param env Pointer to the CPLEX environment.
 @param cbdata Pointer to pass to functions that obtain callback-specific information.
 @param wherefrom An integer value reporting where in the optimization this function was called.
 @param cbhandle Pointer to private user data.
 @param objval_p new objective function value.
 @param x new heuristic solution.
 @param checkfeas_p if 1 CPLEX checks if the solution is feasible, otherwise 0.
 @param useraction_p Pointer to an integer specifying the action for CPLEX to take at the completion of the user callback.
 @return The routine returns 0 (zero) if successful and nonzero if an error occurs.
 */
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
        
        //printf("Obj value: %lf, Thread: %d, Node: %d\n", objval, mythread, node_count);
        
        inst -> flag[mythread] = 0;
        *objval_p = objval;
        *checkfeas_p = 1;
        *useraction_p = CPX_CALLBACK_SET;
    }
    
    return 0;
}

// ********************************************************************************* //

// ******************************** GENERIC CALLBACK ******************************* //

/**
 Generic callback implementation.

 @param context Pointer to a callback context.
 @param contextid context identification.
 @param user Pointer to private user data.
 @return The routine returns 0 (zero) if successful and nonzero if an error occurs.
 */
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
                
                //printf("Obj value: %lf, Thread: %d, Node: %d\n", objval, mythread, node_count);
                if (CPXcallbackpostheursoln(context, inst -> ncols, ind, x, objval, CPXCALLBACKSOLUTION_CHECKFEAS)) print_error("Error in generic callback: Heuristic");
                
                inst -> flag[mythread] = 0;
                free(x);
                free(ind);
            }
            break;
        }
        case CPX_CALLBACKCONTEXT_RELAXATION: {
            
            // total number of nodes solved
            //CPXINT node_count = 0; if (CPXcallbackgetinfoint(context, CPXCALLBACKINFO_NODECOUNT, &node_count)) print_error("Generic callback - separation: get info node count error");
            
            //CPXINT node_depth = 0; if (CPXcallbackgetinfoint(context, CPXCALLBACKINFO_NODEDEPTH, &node_depth)) print_error("Generic callback - separation: get info node depth error");
            //printf("node depth: %d\n", node_depth);
            
            //if (node_depth <= 10) {
            
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
                
                inputGen in;
                in.context = context;
                in.inst = inst;
                
                if(ncomp > 1)
                {
                    // Add constraints on connected components
                    separationMultiple_gen(inst, ncomp, compscount, comps, context);
                }
                if (ncomp == 1) {
                    if ( CCcut_violated_cuts (inst -> nnodes, nedge, elist, xstar, 2.0 - EPSILON, doit_fn_concorde_gen, (void*)&in)) print_error("error in CCcut_violated_cuts");
                }
                
                free(elist);
                free(comps);
                free(compscount);
            
            //}
            
            break;
        }
    }
    
    free(xstar);
    
    return 0;
}

void addSecInRelaxation(CPXCALLBACKCONTEXTptr context, instance* inst, int ncomp, int* compscount, int* comps) {
    
    int cnt = 0;
    char sense = 'L';
    int *index = (int *) calloc((inst->nnodes * inst->nnodes - 1)/2, sizeof(int));
    double *value = (double *) calloc((inst->nnodes * inst->nnodes - 1)/2, sizeof(double));
    int *start_indexes = (int *) calloc((inst->nnodes * inst->nnodes - 1)/2, sizeof(int));
    int type_purge = CPX_USECUT_FILTER;
    
    for (int i = 0; i < ncomp; i++) {
        for (int j = 0; j < compscount[i]; j++) {
            index[j] = comps[cnt++];
            value[j] = 1.0;
        }
        double rhs = compscount[i] - 1;
        if (CPXcallbackaddusercuts(context, 1, rhs, &rhs, &sense, start_indexes, index, value, &type_purge, start_indexes)) print_error("Generic callback - separation: CPXcallbackaddusercuts error in function addSecInRelaxation");
    }
}

/**
 Add SEC constraints to the TSP problem.

 @param inst instance of the struct "instance" for TSP problem.
 @param xstar TSP solution in CPLEX format.
 @param context CPLEX context.
 @return The routine returns 0 (zero) if successful and nonzero if an error occurs.
 */
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
        //CPXINT node_count = 0; if (CPXcallbackgetinfoint(context, CPXCALLBACKINFO_NODECOUNT, &node_count)) print_error("Generic callback - separation: get info node count error");
        
        //CPXINT node_depth = 0; if (CPXcallbackgetinfoint(context, CPXCALLBACKINFO_NODEUID, &node_depth)) print_error("Generic callback - separation: get info node depth error");
        
        // provide to MIP a TSP integer solution only in nodes within depth nine of the branching tree
        // or every 10 nodes
        //if (node_count <= 1024) {
            
            // save the complete graph from components found by CPLEX
            // each thread has a specific solution
            // flag[i] = 1 -> thread i has a solution available
            int mythread = -1; if (CPXcallbackgetinfoint(context, CPXCALLBACKINFO_THREADID, &mythread)) print_error("Generic callback - separation: get info thread error");
            complete_cycle(inst, succ, comp, &ncomp);
            for (int i = 0; i < inst -> nnodes; i++) inst -> sol_thread[mythread][i] = succ[i];
            inst -> flag[mythread] = 1;
        //}
    }
    
    free(succ);
    free(comp);
    
    return (ncomp == 1 ? 0 : ncomp);
}

/**
 Add SEC in continuous relaxation solutions.

 @param cutval value of size of the cut.
 @param cutcount number of nodes in the cut.
 @param cut nodes belonging to the cut.
 @param inParam Pointer to private user data.
 @return The routine returns 0 (zero) if successful and nonzero if an error occurs.
 */
/*
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
*/

/**
 Add SEC in continuous relaxation solutions.
 
 @param cutval value of size of the cut.
 @param cutcount number of nodes in the cut.
 @param cut nodes belonging to the cut.
 @param inParam Pointer to private user data.
 @return The routine returns 0 (zero) if successful and nonzero if an error occurs.
 */
int doit_fn_concorde_gen(double cutval, int cutcount, int *cut , void *inParam) {
    //printf("cutval: %lf \ncutcount: %d \n", cutval, cutcount);
    int i,j;
    
    inputGen* in = (inputGen*) inParam;
    
    int num_x_var = (in -> inst -> nnodes - 1) * in -> inst -> nnodes / 2;     // number of x variables
    int *start_indexes = (int *) calloc(cutcount, sizeof(int));
    int nzcnt = 0;                    // number of non-zero variables in the constraint
    char sense = 'L';
    double* rmatval = (double*) calloc(num_x_var, sizeof(double));        // coefficients of the non-zero variables
    int* rmatind = (int*) calloc(num_x_var, sizeof(int));           // position of the variables to set (in terms of columns)
    double rhs = cutcount - 1;
    for(int i = 0; i < num_x_var; i++)
    {
        rmatval[i] = 0.0;
    }
    
    for(i=0; i < cutcount; ++i)
    {
        for(j=0; j < cutcount; ++j)
        {
            if(cut[i] >= cut[j])
            {
                continue;
            }
            else
            {
                rmatind[nzcnt] = xpos(cut[i], cut[j], in->inst);
                rmatval[nzcnt] = 1.0;
                nzcnt++;
            }
        }
    }
    
    int type_purge = CPX_USECUT_FORCE;
    CPXcallbackaddusercuts(in->context, 1, nzcnt, &rhs, &sense, start_indexes, rmatind, rmatval, &type_purge, start_indexes);
    if(VERBOSE > 50) {
        printf("Cut added with one connected component\n");
    }
    
    free(rmatval);
    free(rmatind);
    free(start_indexes);
    return 0;
}

int separationMultiple_gen(instance *inst, int ncomp, int *compscount, int *comps, CPXCALLBACKCONTEXTptr context)
{
    int num_x_var = (inst -> nnodes-1) * inst -> nnodes / 2;     // number of x variables
    int offset = 0;
    int added = 0;
    int *start_indexes = (int *) calloc((inst -> nnodes * inst -> nnodes - 1)/2, sizeof(int));
    int type_purge = CPX_USECUT_FORCE;
    
    for(int i = 0; i < ncomp; i++)
    {
        int nzcnt = 0;                    // number of non-zero variables in the constraint
        char sense = 'L';
        double* rmatval = (double*) calloc(num_x_var, sizeof(double));        // coefficients of the non-zero variables
        int* rmatind = (int*) calloc(num_x_var, sizeof(int));           // position of the variables to set (in terms of columns)
        
        for(int i = 0; i < num_x_var; i++)
        {
            rmatval[i] = 0;
        }
        
        int comp_vertexes = compscount[i];
        double rhs = comp_vertexes - 1.0;
        
        int j,k;
        for(j = offset; j < offset + comp_vertexes; j++)
        {
            for(k = offset; k < offset + comp_vertexes; k++)
            {
                if((comps[j] == comps[k]) || (comps[j] > comps[k]))
                {
                    continue;
                }
                else
                {
                    rmatind[nzcnt] = xpos(comps[j], comps[k], inst);
                    rmatval[nzcnt] = 1.0;
                    nzcnt++;
                }
            }
            
        }
        
        offset += comp_vertexes;
        if (CPXcallbackaddusercuts(context, 1, nzcnt, &rhs, &sense, start_indexes, rmatind, rmatval, &type_purge, start_indexes)) print_error("Generic callback - separation: CPXcallbackaddusercuts error in function addSecInRelaxation");
        added++;
        
        free(rmatind);
        free(rmatval);
    }
    
    if(VERBOSE > 50) {
        printf("Cut added with multiple connected component\n");
    }
    
    
    free(start_indexes);
    return ncomp;
}

// ********************************************************************************* //

// ****************************** MERGE COMP ALGORITHM ***************************** //

/**
 Build a TSP tour from the distinct components found by a CPLEX solution.
 
 @brief Merge components based on distances between nodes of different components,
 until the solution has one component, that is, a feasible solution.

 @param inst instance of the struct "instance" for TSP problem.
 @param succ TSP solution as successors.
 @param comp number of the component associated to each node in the solution.
 @param ncomp number of components in the solution.
 */
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
        
        for (int i = 0; i < inst -> nnodes - 1; i++) {
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
