//
//  tsp.c
//  cplex
//
//  Created by Matteo Moratello on 13/03/2020.
//  Copyright © 2020 Matteo Moratello. All rights reserved.
//
#include "tsp.h"

void build_model(instance *inst, CPXENVptr env, CPXLPptr lp);
void print_error(const char *err);
void print_solution(instance *inst, double *xsol);
int problem_dimension(instance *inst) { return (inst->nnodes*(inst->nnodes-1))/2; }

// position
int xpos(int i, int j, instance *inst)
{
    if ( i == j ) print_error(" i == j in xpos" );
    if ( i > j ) return xpos(j,i,inst);
    int pos = i * inst->nnodes + j - (( i + 1 ) * ( i + 2 )) / 2;
    return pos;
}

// distance between two nodes
double dist(int i, int j, instance *inst)
{
    double dx = inst->xcoord[i] - inst->xcoord[j];
    double dy = inst->ycoord[i] - inst->ycoord[j];
    int dis = sqrt(dx*dx+dy*dy) + 0.499999999; // nearest integer
    return dis+0.0;
}

// optimizer
int TSPopt(instance *inst) {
    
    // open cplex model
    int error;
    // CPLEX environment
    CPXENVptr env = CPXopenCPLEX(&error);
    // problem object
    CPXLPptr lp = CPXcreateprob(env, &error, "TSP");
    // Cplex's parameter setting
    build_model(inst, env, lp);
    
    // Find solution
    CPXmipopt(env,lp);
    
    double objval;
    // value of objective function
    CPXgetobjval (env, lp, &objval);
    printf("Obj value: %f\n" , objval);
    
    // number of variables of the problem
    int dim = problem_dimension(inst);
    printf("Number of variables: %i\n",problem_dimension(inst));
    
    // final value of variables
    double *bigqstar = (double *) calloc(dim, sizeof(double));
    if (CPXgetx(env, lp, bigqstar, 0, dim-1)) print_error("Error in CPXgetx");
    if ( VERBOSE >= -100 ) for ( int j = 0; j < dim; j++ ) printf(" ... qstar[%3d] = %10.2lf \n", j+1, bigqstar[j]);
    
    // show the solution found
    print_solution(inst, bigqstar);
    
    // free and close cplex model
    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);
    
    return 0;
}

void print_solution(instance *inst, double *xsol) {
    
    FILE *file;
    int count = 0;
    file = fopen("x_variables.txt", "wt");
    
    // save and print edges that correnspond to variable x(*,*) = 1
    for(int i=0; i< inst->nnodes; i++) {
        for(int j=i+1; j< inst->nnodes; j++) {
            if(xsol[count]) {
                printf("Edge x(%i %i)\n", i+1, j+1);
            
                fprintf(file, "%f %f %i\n", inst->xcoord[i], inst->ycoord[i], i+1);
                fprintf(file, "%f %f %i\n\n", inst->xcoord[j], inst->ycoord[j], j+1);
            }
            count ++;
        }
    }
    
    fclose(file);
    
    // set the graph
    file = fopen("solution.txt", "wt");
    
    fprintf(file, "set xrange [-1000:9000] \n set yrange [-1000:6500]\n ");
    //fprintf(file, "set xrange [-300:500] \n set yrange [-300:500]\n ");
    fprintf(file, "set title 'TSP'\n set xlabel 'x'\n set ylabel 'y' \n set style line 5 lt rgb \"#009ad1\" lw 2 pt 6 ps 1 \n");
    fprintf(file, "plot \"x_variables.txt\" with  linespoint ls 5, '' using 1:2:3 with labels tc default font \"arial,12\" offset 0, 1 \n");
    fprintf(file, "pause -1\n");
    fclose(file);
    
    system("/usr/local/Cellar/gnuplot/5.2.8/bin/gnuplot solution.txt");
    
}

// build the model
void build_model(instance *inst, CPXENVptr env, CPXLPptr lp) {
    
    //double zero = 0.0;
    // type: Binary
    char binary = 'B';
    
    char **cname = (char **) calloc(1, sizeof(char*));
    cname[0] = (char *) calloc(100, sizeof(char));
    
    // add binary var.s x(i,j) for i < j
    for ( int i = 0; i < inst->nnodes; i++ )
    {
        for ( int j = i+1; j < inst->nnodes; j++ )
        {
            sprintf(cname[0], "x(%d,%d)", i+1,j+1);
            double obj = dist(i,j,inst); // cost == distance
            double lb = 0.0;
            double ub = 1.0;
            if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname) ) print_error(" wrong CPXnewcols on x var.s");
            if ( CPXgetnumcols(env,lp)-1 != xpos(i,j, inst) ) print_error(" wrong position for x var.s");
        }
    }
    
    // add the degree constraints
    for ( int h = 0; h < inst->nnodes; h++ )         // degree constraints
    {
        int lastrow = CPXgetnumrows(env,lp);
        double rhs = 2.0;
        char sense = 'E';                            // 'E' for equality constraint
        sprintf(cname[0], "degree(%d)", h+1);
        if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) print_error(" wrong CPXnewrows [degree]");
        for ( int i = 0; i < inst->nnodes; i++ )
        {
            if ( i == h ) continue;
            if ( CPXchgcoef(env, lp, lastrow, xpos(i,h, inst), 1.0) ) print_error(" wrong CPXchgcoef [degree]");
        }
    }
    
    // save the model
    if ( VERBOSE >= -100 ) CPXwriteprob(env, lp, "model.lp", NULL);
    
    free(cname[0]);
    free(cname);
    
}
