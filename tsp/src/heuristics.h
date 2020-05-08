//
//  heuristics.h
//  cplex
//
//  Created by Matteo Moratello on 04/05/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//

#ifndef heuristics_h
#define heuristics_h

#include <stdio.h>
#include "tsp.h"

typedef struct tPoint {
    double x, y;
} Point;

int convexHull(Point points[], int n, Point* ch, int* size);

void NearNeigh(instance *inst, double *xstar);
void grasp(instance *inst, double *xstar);
void insertion(instance *inst, double *xstar);
void insertion_ch(instance *inst, double *xstar);

#endif /* heuristics_h */
