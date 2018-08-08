#ifndef _REXstar_JGG_parallel_h
#define _REXstar_JGG_parallel_h

#include <time.h>
#include "Shared.h"

typedef struct{
	int my_rank;
	int n_ps;
	int n_gene;
	int n_generation;
	int n_population;
	int n_parent;
	int n_children;
	double t_rexstar;
	int selection_type;
	double vtr;
	long t_limit;
	int n_constraint;
	double Pf;
	int output_intvl;
	char *out_transition;
	char *out_solution;
	char *out_population;
	FitnessFun fitnessfun;
	DecodingFun decodingfun;
} RCGAParam;

void RCGA(RCGAParam *Param, Chromosome **Population, Chromosome *best, time_t *t0);
void JGG(RCGAParam *Param, Chromosome **Population);
void REXstar(RCGAParam *Param, Chromosome **p, Chromosome **c);
void RCGAInitial(int *argc, char ***argv, int seed, int n_gene, int n_generation, int n_population,
	int n_parent, int n_children, double t_rexstar, int selection_type, double vtr, long t_limit,
	int n_constraint, double Pf, int output_intvl, char *out_transition, char *out_solution, char *out_population,
	FitnessFun fitnessfun, DecodingFun decodingfun,
	RCGAParam **Param, Chromosome ***Population, Chromosome **best, time_t **t0);
void RCGADeInitial(RCGAParam *Param, Chromosome **Population, Chromosome *best, time_t *t0);
void InitPopulation(RCGAParam *Param, Chromosome **Population);

void requestFitnessCalc(RCGAParam *Param, Chromosome **Group, int n_group);
void calcFitness4Master(RCGAParam *Param);
void dismissWorkers(RCGAParam *Param);

#endif
