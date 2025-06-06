/*
Parallel version of REXstar/JGG with stochastic ranking.

libRCGA is distributed under GNU General Public License v3.0.
For academic usage, libRCGA is free.
For other usages, please contact the author(s).


In this program, the below problem will be solved as an example.

--------------------- Example Problem ---------------------
Minimize:
    f(x) = x[1] * x[1] + x[2] * x[2] + ... + x[10] * x[10]
	
Subject to:
    g1(x) = x[1] * x[2] + 1 <= 0
	g2(x) = x[1] + x[2] + 1 <= 0
	-5.12 <= x[i] <= 5.12 for all i
-----------------------------------------------------------

This source code includes the below three functions.
1. main           : RCGA-related parameters are set, and RCGAInitial, RCGA, and RCGADeInitial are called in this order.
2. example        : Problem to be solved is defined.
3. example_decode : Lower and upper bounds of decision variables and how to convert genes to decision variables are defined.
To solve your own problem, you need to modify these three functions.

Compilation:
mpicc example_REXstar_JGG_parallel.c ../lib/librexstarjgg_parallel.a -I../include -o example_REXstar_JGG_parallel_executable -lm

Execution (two processor cores are used):
mpirun -np 2 ./example_REXstar_JGG_parallel_executable

By default, three output files are generated by the execution.
*_*_transition.dat : Best individual during the search. This file is produced every generation.
*_*_solution.dat   : Best individual during the search. This file is produced at the end of the search.
*_*_population.dat : Population at the final generation. This file is produced at the end of the search.

Please note that the example problem is very simple, so parallelization does not speed up the program.
Parallelization is useful when evaluation of f and g's takes a long computational time.
*/

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "REXstar_JGG_parallel.h"
#include "Shared.h"

void example(double *x, int n_gene, int n_constraint, double *f, double *g);
void example_decode(double *gene, int n_gene, double *x);

int main(int argc, char **argv)
{
	RCGAParam *Param;
	Chromosome **Population, *best;
	time_t *t0;
	int n_gene, n_generation, n_population, n_parent, n_children,
		selection_type, n_constraint, output_intvl;
	double t_rexstar, vtr, Pf;
	char *out_transition, *out_solution, *out_population;
	int seed;
	long t_limit;
	FitnessFun fitnessfun;
	DecodingFun decodingfun;
	
	seed = 1;     // Seed of random number generator. This must be positive.
	n_gene = 10;  // Number of decision variables.
	n_generation = 1000; // Maximum number of generations.
	n_population = 50;   // Number of individuals in the population.
	n_parent = n_gene+1; // Number of parents for REXstar. n_gene + 1 is recommended.
	n_children   = 50;   // Number of children.
	t_rexstar    = 6.0;  // t (step size) for REXstar. 2.5 - 15 is recommended.
	selection_type = 0;  // 0: Best individuals are selected from children (Kobayashi 2009), 1: From Family (Kimura at al., 2015)
	vtr = -1e+10; // Objective function value to be reached. Solution is an individual with the objective function (f) <= vtr and the penalty function (phi) = 0.
	t_limit = 60*60;  // Time limit (sec).
	n_constraint = 2; // Number of constraints. If you have any constraint, this must be zero.
	Pf = 0.45;        // Pf for stochastic ranking sort. 0 <= Pf <= 1. Pf matters only if n_constraint >= 1. Pf = 0 indicates naive ranking.
	output_intvl = 1; // Standard output is made every output_intvl generations. If output_intvl < 1, the standard output is not produced.
	out_transition = "REXstar_JGG_parallel_transition.dat"; // File name of the output file which shows the best individual during the search. This file is produced every output_intvl generations. If output_intvl < 1, this output file is not produced.
	out_solution   = "REXstar_JGG_parallel_solution.dat";   // File name of the output file which shows the best individual during the search. This file is produced at the end of the search.
	out_population = "REXstar_JGG_parallel_population.dat"; // File name of the output file which shows the population at the final generation. This file is produced at the end of the search.
	fitnessfun  = &example;       // User-defined function that defines the problem to be solved.
	decodingfun = &example_decode;// User-defined function that defines the lower and upper bounds of decision variables, and how to convert genes to decision variables

	// Initialization.
	RCGAInitial(&argc,&argv,seed,n_gene,n_generation,n_population,
		n_parent,n_children,t_rexstar,selection_type,vtr,t_limit,
		n_constraint,Pf,output_intvl,out_transition,out_solution,out_population,
		fitnessfun, decodingfun,
		&Param,&Population,&best,&t0);

	// Run RCGA. RCGA is terminated when a solution is found, t_limit is elapsed, or n_generation is reached.
	RCGA(Param,Population,best,t0);

	// De-initialization.
	RCGADeInitial(Param,Population,best,t0);

	return(0);
}


void example(double *x, int n_gene, int n_constraint, double *f, double *g)
{
	int i;

	// Objective function.
	// f(x) = x[1] * x[1] + x[2] * x[2] + ... + x[10] * x[10]
	(*f) = 0;
	for(i=1; i<=n_gene; i++) (*f) += x[i] * x[i];

	// Constraints. If you do not have any constraint, do not touch g.
	// g1(x) = x[1] * x[2] + 1 <= 0
	// g2(x) = x[1] + x[2] + 1 <= 0
	g[1] = x[1] * x[2] + 1;
	g[2] = x[1] + x[2] + 1;

	return;
}


void example_decode(double *gene, int n_gene, double *x)
{
	double lb[n_gene+1], ub[n_gene+1];
	int i;

	// Lower and upper bounds.
	// -5.12 <= x[i] <= 5.12 for all i
	for(i=1; i<=n_gene ;i++){
		lb[i] = -5.12;
		ub[i] =  5.12;
	}

	// Decoding gene[i] to x[i].
	// Note that gene[i] takes a value from zero to one.
	for(i=1; i<=n_gene; i++)
		x[i] = gene[i] * ( ub[i] - lb[i] ) + lb[i];

	return;
}
