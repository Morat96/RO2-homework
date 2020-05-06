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

Point* convexHull(Point points[], int n, Point* ch, int* size);

void NearNeigh(instance *inst);
void grasp(instance *inst);
void insertion(instance *inst);
void insertion_ch(instance *inst);

#endif /* heuristics_h */
