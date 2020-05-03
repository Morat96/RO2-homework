//
//  utilities.h
//  cplex
//
//  Created by Matteo Moratello on 13/04/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//

#ifndef utilities_h
#define utilities_h

#include <stdio.h>
#include "tsp.h"

// read and parse input
void read_input(instance *inst);
void parse_command_line(int argc, char** argv, instance *inst);

// print TSP solution
void print_solution(instance *inst, int *succ);

#endif /* utilities_h */
