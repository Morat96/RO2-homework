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

// **************** LAZY CALLBACK *************** //
// lazy callback: integer LP
static int CPXPUBLIC mylazycallback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p);
// add SEC in integer solutions
int myseparation(instance *inst, double *xstar, CPXCENVptr env, void *cbdata, int wherefrom);
// user callback: relaxation
static int CPXPUBLIC UserCutCallback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p);
// add SEC in continuous relaxation solutions
int doit_fn_concorde(double cutval , int cutcount , int *cut , void *inParam);
// ********************************************** //

// ************** GENERIC CALLBACK ************** //
// either integer and relaxation callback
int my_generic_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *user);
// add SEC in integer solutions
int my_separation(instance *inst, double *xstar, CPXCALLBACKCONTEXTptr context);
// add SEC in continuous relaxation solutions
int doit_fn_concorde_gen(double cutval, int cutcount, int *cut , void *inParam);
// ********************************************** //

#endif /* callback_h */
