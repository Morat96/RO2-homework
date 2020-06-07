//
//  callback.h
//  cplex
//
//  Created by Matteo Moratello on 25/04/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//

#ifndef callback_h
#define callback_h

#include <stdio.h>
#include "tsp.h"
#include <cut.h>

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
                             int *useraction_p);

/**
 Add SEC constraints to the problem.
 
 @param inst instance of the struct "instance" for TSP problem.
 @param xstar TSP solution in CPLEX format.
 @param env Pointer to the CPLEX environment.
 @param cbdata Pointer to pass to functions that obtain callback-specific information.
 @param wherefrom An integer value reporting where in the optimization this function was called.
 @return The routine returns 0 (zero) if successful and nonzero if an error occurs.
 */
int myseparation(instance *inst, double *xstar, CPXCENVptr env, void *cbdata, int wherefrom);

/**
 User callback: produce a cut for the LP TSP problem in fractional solutions.
 
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
                              int *useraction_p);

/**
 Add SEC in continuous relaxation solutions.
 
 @param cutval value of size of the cut.
 @param cutcount number of nodes in the cut.
 @param cut nodes belonging to the cut.
 @param inParam Pointer to private user data.
 @return The routine returns 0 (zero) if successful and nonzero if an error occurs.
 */
int doit_fn_concorde(double cutval , int cutcount , int *cut , void *inParam);

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
                           int        *useraction_p);

/**
 Generic callback implementation.
 
 @param context Pointer to a callback context.
 @param contextid context identification.
 @param user Pointer to private user data.
 @return The routine returns 0 (zero) if successful and nonzero if an error occurs.
 */
int my_generic_callback(CPXCALLBACKCONTEXTptr context,
                        CPXLONG contextid,
                        void *user);

/**
 Add SEC constraints to the TSP problem.
 
 @param inst instance of the struct "instance" for TSP problem.
 @param xstar TSP solution in CPLEX format.
 @param context CPLEX context.
 @return The routine returns 0 (zero) if successful and nonzero if an error occurs.
 */
int my_separation(instance *inst, double *xstar, CPXCALLBACKCONTEXTptr context);

/**
 Add SEC in continuous relaxation solutions.
 
 @param cutval value of size of the cut.
 @param cutcount number of nodes in the cut.
 @param cut nodes belonging to the cut.
 @param inParam Pointer to private user data.
 @return The routine returns 0 (zero) if successful and nonzero if an error occurs.
 */
int doit_fn_concorde_gen(double cutval, int cutcount, int *cut , void *inParam);

/**
 Build a TSP tour from the distinct components found by a CPLEX solution.
 
 @brief Merge components based on distances between nodes of different components,
 until the solution has one component, that is, a feasible solution.
 
 @param inst instance of the struct "instance" for TSP problem.
 @param succ TSP solution as successors.
 @param comp number of the components of each node in the solution.
 @param ncomp number of components in the solution.
 */
void complete_cycle(instance *inst, int *succ, int *comp, int *ncomp);

//////////////////////////// utils ////////////////////////////

/**
 Print an error and exit from the program.
 
 @param err error to show.
 */
void print_error(const char *err);

/**
 Return the index of CPLEX solution array from two subsequent nodes..
 
 @param i first node.
 @param j second node.
 @param inst instance of the struct "instance" for TSP problem.
 @return array position.
 */
int xpos(int i, int j, instance *inst);

/**
 @brief Compute the distance between two nodes.
 
 @param i first node.
 @param j second node.
 @param inst instance of the struct "instance" for TSP problem.
 @return distance between two nodes.
 */
double dist(int i, int j, instance *inst);

/**
 Build succ() and comp() wrt xstar().
 
 @param xstar CPLEX solution.
 @param inst instance of the struct "instance" for TSP problem.
 @param succ TSP solution as successors.
 @param comp number of the component associated to each node in the solution.
 @param ncomp number of components in the solution.
 */
void build_sol(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp);

#endif /* callback_h */
