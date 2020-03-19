//
//  main.c
//  cplex
//
//  Created by Matteo Moratello on 12/03/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//

#include <cplex.h>
#include <stdio.h>
#include "tsp.h"

void read_input(instance *inst);
void parse_command_line(int argc, char** argv, instance *inst);
void print_error(const char *err);

void print_error(const char *err) { printf("\n\n ERROR: %s \n\n", err); fflush(NULL); exit(1); }

void free_instance(instance *inst)
{
    free(inst->demand);
    free(inst->xcoord);
    free(inst->ycoord);
}

float f(float x) {
    return(sin(x)*15 + x*x);
}

void print_vertices(instance *inst);

int main(int argc, char **argv) {
    
    if ( argc < 2 ) { printf("Usage: %s -help for help\n", argv[0]); exit(1); }
    if ( VERBOSE >= 2 ) { for (int a = 0; a < argc; a++) printf("%s ", argv[a]); printf("\n"); }
    
    instance inst;
    
    parse_command_line(argc, argv, &inst);
    read_input(&inst);
    print_vertices(&inst);
    free_instance(&inst);
    
    return 0;
}

void print_vertices(instance *inst) {
    
    FILE *file;
    file = fopen("vertices.txt", "wt");
    for(int i=0; i< inst->nnodes; i++) fprintf(file, "%f %f %i\n", inst->xcoord[i], inst->ycoord[i], i+1);
    fclose(file);
    
    file = fopen("command.txt", "wt");
    
    fprintf(file, "set xrange [-1000:9000] \n set yrange [-1000:6500]\n ");
    fprintf(file, "set title 'TSP'\n set xlabel 'x'\n set ylabel 'y' \n set style line 5 lt rgb \"#009ad1\" lw 2 pt 6 ps 1 \n");
    fprintf(file, "plot \"vertices.txt\" with  linespoint ls 5, '' using 1:2:3 with labels tc default font \"arial,12\" offset 0, 1 \n");
    fprintf(file, "pause -1\n");
    fclose(file);
    
    system("gnuplot command.txt");
    
}

// simplified Parser
void read_input(instance *inst) {
    
    FILE *fin = fopen(inst->input_file, "r");
    if ( fin == NULL ) print_error(" input file not found!");
    
    inst->nnodes = -1;
    
    char line[180];
    char *par_name;
    char *token1;
    char *token2;
    
    int active_section = 0; // =1 NODE_COORD_SECTION
    
    int do_print = ( VERBOSE >= 500 );
    
    while ( fgets(line, sizeof(line), fin) != NULL )
    {
        if ( VERBOSE >= 2000 ) { printf("%s",line); fflush(NULL); }
        if ( strlen(line) <= 1 ) continue; // skip empty lines
        par_name = strtok(line, " :");
        if ( VERBOSE >= 3000 ) { printf("parameter \"%s\" ",par_name); fflush(NULL); }
        
        if ( strncmp(par_name, "NAME", 4) == 0 )
        {
            active_section = 0;
            continue;
        }
        
        if ( strncmp(par_name, "COMMENT", 7) == 0 )
        {
            active_section = 0;
            token1 = strtok(NULL, "");
            if ( VERBOSE >= 10 ) printf(" ... solving instance %s with model %d\n\n", token1, inst->model_type);
            continue;
        }
        
        if ( strncmp(par_name, "TYPE", 4) == 0 )
        {
            token1 = strtok(NULL, " :");
            //if ( strncmp(token1, "CVRP",4) != 0 ) print_error(" format error:  only TYPE == CVRP implemented so far!!!!!!");
            active_section = 0;
            continue;
        }
        
        if ( strncmp(par_name, "DIMENSION", 9) == 0 )
        {
            if ( inst->nnodes >= 0 ) print_error(" repeated DIMENSION section in input file");
            token1 = strtok(NULL, " :");
            inst->nnodes = atoi(token1);
            if ( do_print ) printf(" ... nnodes %d\n", inst->nnodes);
            inst->demand = (double *) calloc(inst->nnodes, sizeof(double));
            inst->xcoord = (double *) calloc(inst->nnodes, sizeof(double));
            inst->ycoord = (double *) calloc(inst->nnodes, sizeof(double));
            active_section = 0;
            continue;
        }
        
        if ( strncmp(par_name, "EDGE_WEIGHT_TYPE", 16) == 0 )
        {
            token1 = strtok(NULL, " :");
            //if ( strncmp(token1, "EUC_2D", 6) != 0 ) print_error(" format error:  only EDGE_WEIGHT_TYPE == EUC_2D implemented so far!!!!!!");
            active_section = 0;
            continue;
        }
        
        if ( strncmp(par_name, "NODE_COORD_SECTION", 18) == 0 )
        {
            if ( inst->nnodes <= 0 ) print_error(" ... DIMENSION section should appear before NODE_COORD_SECTION section");
            active_section = 1;
            continue;
        }
        
        if ( strncmp(par_name, "EOF", 3) == 0 )
        {
            break;
        }
        
        if ( active_section == 1 ) // within NODE_COORD_SECTION
        {
            int i = atoi(par_name) - 1;
            if ( i < 0 || i >= inst->nnodes ) print_error(" ... unknown node in NODE_COORD_SECTION section");
            token1 = strtok(NULL, " :,");
            token2 = strtok(NULL, " :,");
            inst->xcoord[i] = atof(token1);
            inst->ycoord[i] = atof(token2);
            if ( do_print ) printf(" ... node %4d at coordinates ( %15.7lf , %15.7lf )\n", i+1, inst->xcoord[i], inst->ycoord[i]);
            continue;
        }
        
        printf(" final active section %d\n", active_section);
        print_error(" ... wrong format for the current simplified parser!!!!!!!!!");
        
    }
    fclose(fin);
    
    // print information about nodes
    if(VERBOSE >= 50) {
        printf("RECAP\nTotal nodes: %d\n", inst -> nnodes);
        printf(" ------------------------------------\n");
        printf("| %4s |  %6s      |  %6s      |\n","node","x","y");
        printf("|------------------------------------|\n");
        for(int i=0; i< inst->nnodes; i++) printf("|%4d  |%13.6lf |%13.6lf |\n", i+1, inst->xcoord[i], inst->ycoord[i]);
        printf(" ------------------------------------\n");
    }
}

// parse command line
void parse_command_line(int argc, char** argv, instance *inst)
{
    if ( VERBOSE >= 100 ) printf(" running %s with %d parameters \n", argv[0], argc-1);
    
    // default
    inst->model_type = 0;
    strcpy(inst->input_file, "NULL");
    inst->timelimit = CPX_INFBOUND;
    
    int help = 0; if ( argc < 1 ) help = 1;
    for ( int i = 1; i < argc; i++ )
    {
        if ( strcmp(argv[i],"-file") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; }       // input file
        if ( strcmp(argv[i],"-time_limit") == 0 ) { inst->timelimit = atof(argv[++i]); continue; }  // total time limit
        if ( strcmp(argv[i],"-model_type") == 0 ) { inst->model_type = atoi(argv[++i]); continue; } // model type
        help = 1;
    }
    
    if ( help || (VERBOSE >= 10) )        // print current parameters
    {
        printf("\n\navailable parameters (vers. 13-mar-2020) --------------------------------------------------\n");
        printf("-file %s\n", inst->input_file);
        printf("-time_limit %lf\n", inst->timelimit);
        printf("-model_type %d\n", inst->model_type);
        printf("-------------------------------------------------------------------------------------------\n\n");
    }
    
    if ( help ) exit(1);
    
}

