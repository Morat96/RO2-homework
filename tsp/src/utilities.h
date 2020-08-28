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

/**
 Input parser.
 
 @param inst instance of the struct "instance" for TSP problem.
 */
void read_input(instance *inst);

void read_input1(char **filename ,instance *inst);

/**
 Parse the command line.
 
 @param argc number of arguments.
 @param argv array of string with arguments.
 @param inst instance of the struct "instance" for TSP problem.
 */
void parse_command_line(int argc, char** argv, instance *inst);

/**
 Print the tour.
 
 @param inst instance of the struct "instance" for TSP problem.
 @param succ tour solution described with successors.
 */
void print_solution(instance *inst, int *succ);

/**
 Print a light version of the tour.
 
 @param inst instance of the struct "instance" for TSP problem.
 @param succ tour solution described with successors.
 */
void print_solution_light(instance *inst, int *succ);

/**
 Print an error and exit from the program.
 
 @param err error to show.
 */
void print_error(const char *err);

#endif /* utilities_h */
