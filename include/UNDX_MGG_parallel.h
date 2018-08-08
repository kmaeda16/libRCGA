#ifndef _UNDX_MGG_parallel_h
#define _UNDX_MGG_parallel_h

#include <time.h>
#include "Shared.h"

typedef struct{
	int my_rank;
	int n_ps;
	int n_gene;
	int n_generation;
	int n_population;
	int n_children;
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
void MGG(RCGAParam *Param, Chromosome **Population);
void getNewChild(RCGAParam *Param, Chromosome *p1, Chromosome *p2, Chromosome *p3, Chromosome *c);
void UNDX(RCGAParam *Param, Chromosome *p1, Chromosome *p2, Chromosome *p3, Chromosome *c);
void RCGAInitial(int *argc, char ***argv, int seed, int n_gene, int n_generation, int n_population,
	int n_children, double vtr, long t_limit,
	int n_constraint, double Pf, int output_intvl, char *out_transition, char *out_solution, char *out_population,
	FitnessFun fitnessfun, DecodingFun decodingfun,
	RCGAParam **Param, Chromosome ***Population, Chromosome **best, time_t **t0);
void RCGADeInitial(RCGAParam *Param, Chromosome **Population, Chromosome *best, time_t *t0);
void InitPopulation(RCGAParam *Param, Chromosome **Population);

void requestFitnessCalc(RCGAParam *Param, Chromosome **Group, int n_group);
void calcFitness4Master(RCGAParam *Param);
void dismissWorkers(RCGAParam *Param);

#endif
