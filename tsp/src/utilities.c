//
//  utilities.c
//  cplex
//
//  Created by Matteo Moratello on 13/04/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//

#include "utilities.h"

/**
 Input parser.

 @param inst instance of the struct "instance" for TSP problem.
 */
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
            if ( strncmp(token1, "TSP", 3) != 0 ) print_error(" format error:  only TYPE == TSP implemented so far!!!!!!");
            active_section = 0;
            continue;
        }
        
        if ( strncmp(par_name, "DIMENSION", 9) == 0 )
        {
            if ( inst->nnodes >= 0 ) print_error(" repeated DIMENSION section in input file");
            token1 = strtok(NULL, " :");
            inst->nnodes = atoi(token1);
            if ( do_print ) printf(" ... nnodes %d\n", inst->nnodes);
            inst->xcoord = (double *) calloc(inst->nnodes, sizeof(double));
            inst->ycoord = (double *) calloc(inst->nnodes, sizeof(double));
            inst -> sol_thread = (int **) calloc(4, sizeof(int*));
            inst -> sol_thread[0] = (int *) calloc(inst -> nnodes, sizeof(int));
            inst -> sol_thread[1] = (int *) calloc(inst -> nnodes, sizeof(int));
            inst -> sol_thread[2] = (int *) calloc(inst -> nnodes, sizeof(int));
            inst -> sol_thread[3] = (int *) calloc(inst -> nnodes, sizeof(int));
            for (int t = 0; t < 4; t++) inst -> flag[t] = 0;
            active_section = 0;
            continue;
        }
        
        if ( strncmp(par_name, "EDGE_WEIGHT_TYPE", 16) == 0 )
        {
            token1 = strtok(NULL, " :");
            if ( strncmp(token1, "EUC_2D", 6) == 0 ) inst -> euc_2d = 1;
            else if ( strncmp(token1, "ATT", 3) == 0 ) inst -> att = 1;
            else if ( strncmp(token1, "GEO", 3) == 0 ) inst -> geo = 1;
            else print_error(" format error:  only EDGE_WEIGHT_TYPE == EUC_2D/ATT/GEO implemented so far!!!!!!");
            active_section = 0;
            continue;
        }
        
        if ( strncmp(par_name, "EDGE_WEIGHT_FORMAT", 18) == 0 )
        {
            active_section = 0;
            continue;
        }
        
        if ( strncmp(par_name, "DISPLAY_DATA_TYPE", 17) == 0 )
        {
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

/**
 Input parser.
 
 @param inst instance of the struct "instance" for TSP problem.
 */
void read_input1(char **filename ,instance *inst) {
    
    FILE *fin = fopen(*filename, "r");
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
            if ( strncmp(token1, "TSP", 3) != 0 ) print_error(" format error:  only TYPE == TSP implemented so far!!!!!!");
            active_section = 0;
            continue;
        }
        
        if ( strncmp(par_name, "DIMENSION", 9) == 0 )
        {
            if ( inst->nnodes >= 0 ) print_error(" repeated DIMENSION section in input file");
            token1 = strtok(NULL, " :");
            inst->nnodes = atoi(token1);
            if ( do_print ) printf(" ... nnodes %d\n", inst->nnodes);
            inst->xcoord = (double *) calloc(inst->nnodes, sizeof(double));
            inst->ycoord = (double *) calloc(inst->nnodes, sizeof(double));
            inst -> sol_thread = (int **) calloc(4, sizeof(int*));
            inst -> sol_thread[0] = (int *) calloc(inst -> nnodes, sizeof(int));
            inst -> sol_thread[1] = (int *) calloc(inst -> nnodes, sizeof(int));
            inst -> sol_thread[2] = (int *) calloc(inst -> nnodes, sizeof(int));
            inst -> sol_thread[3] = (int *) calloc(inst -> nnodes, sizeof(int));
            for (int t = 0; t < 4; t++) inst -> flag[t] = 0;
            active_section = 0;
            continue;
        }
        
        if ( strncmp(par_name, "EDGE_WEIGHT_TYPE", 16) == 0 )
        {
            token1 = strtok(NULL, " :");
            inst -> euc_2d = 0;
            inst -> att = 0;
            inst -> geo = 0;
            if ( strncmp(token1, "EUC_2D", 6) == 0 ) inst -> euc_2d = 1;
            else if ( strncmp(token1, "ATT", 3) == 0 ) inst -> att = 1;
            else if ( strncmp(token1, "GEO", 3) == 0 ) inst -> geo = 1;
            else print_error(" format error:  only EDGE_WEIGHT_TYPE == EUC_2D/ATT/GEO implemented so far!!!!!!");
            active_section = 0;
            continue;
        }
        
        if ( strncmp(par_name, "EDGE_WEIGHT_FORMAT", 18) == 0 )
        {
            active_section = 0;
            continue;
        }
        
        if ( strncmp(par_name, "DISPLAY_DATA_TYPE", 17) == 0 )
        {
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

/**
 Parse the command line.

 @param argc number of arguments.
 @param argv array of string with arguments.
 @param inst instance of the struct "instance" for TSP problem.
 */
void parse_command_line(int argc, char** argv, instance *inst)
{
    if ( VERBOSE >= 100 ) printf(" running %s with %d parameters \n", argv[0], argc-1);
    
    // default
    inst -> model_type = 0;
    strcpy(inst -> input_file, "NULL");
    inst -> timelimit = CPX_INFBOUND;
    inst -> randomseed = 201909284;
    inst -> loop = 0;
    inst -> callback = 0;
    inst -> hardfixing = 0;
    inst -> localbranching = 0;
    inst -> vns = 0;
    inst -> tabu_search = 0;
    inst -> sim_annealing = 0;
    inst -> genetic = 0;
    
    int help = 0; if ( argc < 1 ) help = 1;
    for ( int i = 1; i < argc; i++ )
    {
        if ( strcmp(argv[i],"-file") == 0 ) { strcpy(inst -> input_file,argv[++i]); continue; }               // input file
        if ( strcmp(argv[i],"-time_limit") == 0 ) { inst -> timelimit = atof(argv[++i]); continue; }          // total time limit
        if ( strcmp(argv[i],"-model_type") == 0 ) { inst -> model_type = atoi(argv[++i]); continue; }         // model type
        if ( strcmp(argv[i],"-randomseed") == 0 ) { inst -> randomseed = atoi(argv[++i]); continue; }         // change randomseed
        if ( strcmp(argv[i],"-loop") == 0 ) { inst -> loop = atoi(argv[++i]); continue; }                     // loop method
        if ( strcmp(argv[i],"-callback") == 0 ) { inst -> callback = atoi(argv[++i]); continue; }             // callback method
        if ( strcmp(argv[i],"-hardfixing") == 0 ) { inst -> hardfixing = atoi(argv[++i]); continue; }         // hard fixing method
        if ( strcmp(argv[i],"-localbranching") == 0 ) { inst -> localbranching = atoi(argv[++i]); continue; } // local branching method
        if ( strcmp(argv[i],"-vns") == 0 ) { inst -> vns = atoi(argv[++i]); continue; }                       // VNS
        if ( strcmp(argv[i],"-ts") == 0 ) { inst -> tabu_search = atoi(argv[++i]); continue; }                // Tabù Search
        if ( strcmp(argv[i],"-sa") == 0 ) { inst -> sim_annealing = atoi(argv[++i]); continue; }              // Simulated Annealing
        if ( strcmp(argv[i],"-ga") == 0 ) { inst -> genetic = atoi(argv[++i]); continue; }                    // Genetic algorithm
        help = 1;
    }
    
    if ( help || (VERBOSE >= 10) )        // print current parameters
    {
        printf("\n\navailable parameters (vers. 6-jun-2020) --------------------------------------------------\n");
        printf("-file %s\n", inst -> input_file);
        printf("-time_limit %lf\n", inst -> timelimit);
        printf("-model_type %d\n", inst -> model_type);
        printf("-randomseed %d\n", inst -> randomseed);
        printf("-loop %d\n", inst -> loop);
        printf("-callback %d\n", inst -> callback);
        printf("-hardfixing %d\n", inst -> hardfixing);
        printf("-localbranching %d\n", inst -> localbranching);
        printf("-vns %d\n", inst -> vns);
        printf("-ts %d\n", inst -> tabu_search);
        printf("-sa %d\n", inst -> sim_annealing);
        printf("-ga %d\n", inst -> genetic);
        printf("-------------------------------------------------------------------------------------------\n\n");
    }
    
    if ( help ) exit(1);
}

/**
 Print the tour.

 @param inst instance of the struct "instance" for TSP problem.
 @param succ tour solution described with successors.
 */
void print_solution(instance *inst, int *succ) {
    
    FILE *file;
    
    file = fopen("xsol.txt", "wt");
    
    if(VERBOSE >= 50) printf("The solution includes the following edges:\n");
    // save and print edges that correnspond to variable x(*,*) = 1
    for(int i=0; i < inst->nnodes; i++) {
        if(VERBOSE >= 50) printf("Edge x[%i,%i]\n", i+1, succ[i]+1);
        // save edge (i,j)
        fprintf(file, "%f %f %i\n", inst->xcoord[i], inst->ycoord[i], i+1);
        fprintf(file, "%f %f %i\n\n", inst->xcoord[succ[i]], inst->ycoord[succ[i]], succ[i]+1);
    }
    
    fclose(file);
    
    system("/usr/local/Cellar/gnuplot/5.2.8/bin/gnuplot script.txt");
    
}

/**
 Print a light version of the tour.
 
 @param inst instance of the struct "instance" for TSP problem.
 @param succ tour solution described with successors.
 */
void print_solution_light(instance *inst, int *succ) {
    
    FILE *file;
    
    file = fopen("xsol.txt", "wt");
    
    if(VERBOSE >= 50) printf("The solution includes the following edges:\n");
    // save and print edges that correnspond to variable x(*,*) = 1
    for(int i=0; i < inst->nnodes; i++) {
        if(VERBOSE >= 50) printf("Edge x[%i,%i]\n", i+1, succ[i]+1);
        // save edge (i,j)
        fprintf(file, "%f %f\n", inst->xcoord[i], inst->ycoord[i]);
        fprintf(file, "%f %f\n\n", inst->xcoord[succ[i]], inst->ycoord[succ[i]]);
    }
    
    fclose(file);
    
    system("/usr/local/Cellar/gnuplot/5.2.8/bin/gnuplot script_light.txt");
    
}
