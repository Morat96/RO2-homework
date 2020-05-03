//
//  hardfixing.h
//  cplex
//
//  Created by Matteo Moratello on 26/04/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//

#ifndef hardfixing_h
#define hardfixing_h

#include <stdio.h>
#include "tsp.h"

void hardfixing(instance *inst, CPXENVptr env, CPXLPptr lp);

#endif /* hardfixing_h */
