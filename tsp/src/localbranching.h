//
//  localbranching.h
//  cplex
//
//  Created by Matteo Moratello on 27/04/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//

#ifndef localbranching_h
#define localbranching_h

#include <stdio.h>
#include "tsp.h"

void localbranching(instance *inst, CPXENVptr env, CPXLPptr lp);

#endif /* localbranching_h */
