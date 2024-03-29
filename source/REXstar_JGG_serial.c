#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "REXstar_JGG_serial.h"


void RCGA(RCGAParam *Param, Chromosome **Population, Chromosome *best, time_t *t0)
{
	int i, index;
	int n_population, n_gene, n_constraint, n_generation, output_intvl;
	double vtr;
	long t_limit;
	char *out_transition;
	DecodingFun decodingfun;
	int flg_printed = 0;

	n_population = Param->n_population;
	n_gene = Param->n_gene;
	n_constraint = Param->n_constraint;
	n_generation = Param->n_generation;
	output_intvl = Param->output_intvl;
	vtr = Param->vtr;
	t_limit = Param->t_limit;
	out_transition = Param->out_transition;
	decodingfun = Param->decodingfun;

	i = 1;
	InitPopulation(Param,Population);
	index = findBest(Population,n_population);
	copyChrom(best,Population[index],n_gene,n_constraint);

	if(0 < output_intvl){
		printTransition(t0,i,best->f,best->phi);
		writeTransition(t0,i,out_transition,best,n_gene,n_constraint,decodingfun);
		flg_printed = 1;
	}
	if( (best->phi == 0 && best->f <= vtr) || getElapsedTime(t0) >= t_limit ) goto end;
	
	while(i < n_generation){
		i++;
		flg_printed = 0;
		JGG(Param,Population);
		index = findBest(Population,n_population);
		if( Population[index]->phi < best->phi 
			|| (Population[index]->phi == best->phi && Population[index]->f < best->f) )
			copyChrom(best,Population[index],n_gene,n_constraint);
		if(0 < output_intvl && i % output_intvl == 0){
			printTransition(t0,i,best->f,best->phi);
			writeTransition(t0,i,out_transition,best,n_gene,n_constraint,decodingfun);
			flg_printed = 1;
		}
		if( (best->phi == 0 && best->f <= vtr) || getElapsedTime(t0) >= t_limit ) goto end;
	}

end:
	if(0 < output_intvl && flg_printed == 0){
		printTransition(t0,i,best->f,best->phi);
		writeTransition(t0,i,out_transition,best,n_gene,n_constraint,decodingfun);
	}

	return;
}


void JGG(RCGAParam *Param, Chromosome **Population)
{
	/*
	Just Generation Gap Selection for REXstar
	This function alters population using JGG with REXstar.

	See "Kimura S, Sato M, Okada-Hatakeyama M: An Effective Method for the 
	Inference of Reduced S-system Models of Genetic Networks. Information and 
	Media Technologies 2015, 10(1):166-174.".
	*/

	const int maxitr = 10;
	int *ip;
	Chromosome **p, **c, **f;
	int n_parent, n_children, n_gene, n_constraint, n_population, selection_type;
	int flg;
	int i, j, k;
	double Pf;

	n_parent = Param->n_parent;
	n_children = Param->n_children;
	n_gene = Param->n_gene;
	n_constraint = Param->n_constraint;
	n_population = Param->n_population;
	selection_type = Param->selection_type;
	Pf = Param->Pf;

	ip = allocIVector(n_parent);
	p = allocPopulation(n_parent,n_gene,n_constraint);
	c = allocPopulation(n_children,n_gene,n_constraint);

	// Pick up parents from main population
	for(i=1; i<=n_parent; i++){
		for(j=1; j<=maxitr; j++){
			flg = 0;
			ip[i] = (int)(Rand()*(double)n_population) + 1;
			for(k=1; k<i; k++){
				if(ip[k] == ip[i]) flg = 1;
			}
			if(flg == 0) break;
		}
	
	}
	for(i=1; i<=n_parent; i++) copyChrom(p[i],Population[ip[i]],n_gene,n_constraint);

	REXstar(Param,p,c);

	switch(selection_type){
		case 0:
			// Chosen from children (Kobayashi, 2009)
			for(i=1; i<=n_parent; i++) copyChrom(Population[ip[i]],c[i],n_gene,n_constraint);
			break;
		case 1:
			// Chosen from family (Kimura et al., 2015)
			f = allocPopulation(n_parent+n_children,n_gene,n_constraint);
			for(i=1; i<=n_parent; i++) copyChrom(f[i],p[i],n_gene,n_constraint);
			for(i=1; i<=n_children; i++) copyChrom(f[n_parent+i],c[i],n_gene,n_constraint);
			sort(f,n_parent+n_children,n_gene,n_constraint,Pf);
			for(i=1; i<=n_parent; i++) copyChrom(Population[ip[i]],f[i],n_gene,n_constraint);
			freePopulation(f,n_parent+n_children,n_gene,n_constraint);
			break;
		default:
			handleError("Unexpected selection_type!");
			exit(1);
	}
	sort(Population,n_population,n_gene,n_constraint,Pf);

	freePopulation(p,n_parent,n_gene,n_constraint);
	freePopulation(c,n_children,n_gene,n_constraint);
	freeIVector(ip,n_parent);

	return;
}


void REXstar(RCGAParam *Param, Chromosome **p, Chromosome **c)
{
	/*
	Function of REXstar (Real-coded Ensemble Crossover star)
	This function returns children generated by REXstar.

	See "Kimura S, Sato M, Okada-Hatakeyama M: An Effective Method for the 
	Inference of Reduced S-system Models of Genetic Networks. Information and 
	Media Technologies 2015, 10(1):166-174.".

	Recommended:
	n_population = 2*n_gene ~ 20*n_gene
	n_parent = n_gene + 1
	n_children = 2*n_gene ~ 3*n_gene
	t (Step-size parameter) = 2.5 ~ 15

	Required:
	n_parent <= n_children
	*/

	const int maxitr = 100;
	Chromosome *G, *Gb;
	Chromosome **r;
	double xi;
	int flg_oob, flg_eval;
	int i, j, k, l;
	int n_gene, n_parent, n_children, n_constraint;
	double t;
	double Pf;

	n_gene = Param->n_gene;
	n_parent = Param->n_parent;
	n_children = Param->n_children;
	n_constraint = Param->n_constraint;
	t = Param->t_rexstar;
	Pf = Param->Pf;

	if(n_children < n_parent){
		handleError("n_parent <= n_children must be satisfied!\n");
		exit(1);
	}
	
	// Calculate center of gravity G
	G = allocChrom(n_gene,n_constraint);
	for(j=1; j<=n_gene; j++) G->gene[j] = 0.0;
	G->f = G->phi = DBL_MAX;
	for(i=1; i<=n_parent; i++){
		for(j=1; j<=n_gene; j++) G->gene[j] += p[i]->gene[j];
	}
	for(j=1; j<=n_gene; j++) G->gene[j] /= (double)n_parent;

	// Initialize children
	for(i=1; i<=n_children; i++){
		for(j=1; j<=n_gene; j++) c[i]->gene[j] = DBL_MAX;
		c[i]->f = c[i]->phi = DBL_MAX;
	}

	// Calculate reflection points (which are treated as the children)
	for(i=1; i<=n_parent; i++){
		for(j=1; j<=n_gene; j++) c[i]->gene[j] = 2.0 * G->gene[j] - p[i]->gene[j];
		flg_oob = 0;
		for(j=1; j<=n_gene; j++){
			if(c[i]->gene[j] < 0.0 || 1.0 < c[i]->gene[j]) flg_oob = 1;
		}
		if(flg_oob == 0){
			setFitness(c[i],Param->n_gene,Param->n_constraint,Param->fitnessfun,Param->decodingfun);
		}else{
			c[i]->f = c[i]->phi = DBL_MAX;
		}
	}

	// Store parents and refrection points
	r = allocPopulation(2*n_parent,n_gene,n_constraint);
	for(i=1; i<=n_parent; i++){
		copyChrom(r[i],p[i],n_gene,n_constraint);
		copyChrom(r[n_parent+i],c[i],n_gene,n_constraint);
	}
	sort(r,2*n_parent,n_gene,n_constraint,Pf);
	
	// Calculate center of gravity Gb
	Gb = allocChrom(n_gene,n_constraint);
	for(j=1; j<=n_gene; j++) Gb->gene[j] = 0.0;
	Gb->f = Gb->phi = DBL_MAX;
	for(i=1; i<=n_parent; i++){
		for(j=1; j<=n_gene; j++) Gb->gene[j] += r[i]->gene[j];
	}
	for(j=1; j<=n_gene; j++) Gb->gene[j] /= (double)n_parent;

	// Generate children
	for(i=1; i<=n_children; i++){
		flg_eval = 0;
		// If the generated child is out of range, re-generate a child
		for(l=0; l<maxitr; l++){
			flg_oob = 0;
			for(j=1; j<=n_gene; j++){
				if(c[i]->gene[j] < 0.0 || 1.0 < c[i]->gene[j]) flg_oob = 1;
			}
			if(flg_oob == 0) break;
			// Generation
			for(j=1; j<=n_gene; j++) c[i]->gene[j] = G->gene[j];
			for(j=1; j<=n_gene; j++){
				xi = t * Rand();
				c[i]->gene[j] += xi * (Gb->gene[j] - G->gene[j]);
			}
			for(k=1; k<=n_parent; k++){
				xi = 2.0 * sqrt(3.0/((double)n_parent)) * Rand() - sqrt(3.0/((double)n_parent));
				for(j=1; j<=n_gene; j++) c[i]->gene[j] += xi * (p[k]->gene[j] - G->gene[j]);
			}
			flg_eval = 1;
		}
		// If the generated child is out of range, alter the child
		if(flg_oob == 1){
			for(j=1; j<=n_gene; j++){
				if(1.0 < c[i]->gene[j]){
					c[i]->gene[j] = 1.0;
				}else if(c[i]->gene[j] < 0.0){
					c[i]->gene[j] = 0.0;
				}
			}
			flg_eval = 1;
		}
		if(flg_eval == 1) setFitness(c[i],Param->n_gene,Param->n_constraint,Param->fitnessfun,Param->decodingfun);
	}
	sort(c,n_children,n_gene,n_constraint,Pf);
	
	freeChrom(G,n_gene,n_constraint);
	freeChrom(Gb,n_gene,n_constraint);
	freePopulation(r,2*n_parent,n_gene,n_constraint);

	return;
}


void RCGAInitial(int *argc, char ***argv, int seed, int n_gene, int n_generation, int n_population,
	int n_parent, int n_children, double t_rexstar, int selection_type, double vtr, long t_limit,
	int n_constraint, double Pf, int output_intvl, char *out_transition, char *out_solution, char *out_population,
	FitnessFun fitnessfun, DecodingFun decodingfun,
	RCGAParam **Param, Chromosome ***Population, Chromosome **best, time_t **t0)
{
	srand(seed);

	(*Param) = (RCGAParam *)malloc((size_t)sizeof(RCGAParam));
	(*Param)->n_gene = n_gene;
	(*Param)->n_generation = n_generation;
	(*Param)->n_population = n_population;
	(*Param)->n_parent = n_parent;
	(*Param)->n_children = n_children;
	(*Param)->t_rexstar = t_rexstar;
	(*Param)->selection_type = selection_type;
	(*Param)->vtr = vtr;
	(*Param)->t_limit = t_limit;
	(*Param)->n_constraint = n_constraint;
	(*Param)->Pf = Pf;
	(*Param)->output_intvl = output_intvl;
	(*Param)->out_transition = out_transition;
	(*Param)->out_solution = out_solution;
	(*Param)->out_population = out_population;
	(*Param)->fitnessfun = fitnessfun;
	(*Param)->decodingfun = decodingfun;

	(*Population) = allocPopulation(n_population,n_gene,n_constraint);
	(*best) = allocChrom(n_gene,n_constraint);
	(*t0) = (time_t *)malloc((size_t)sizeof(time_t));
	(**t0) = time(NULL);

	return;
}


void RCGADeInitial(RCGAParam *Param, Chromosome **Population, Chromosome *best, time_t *t0)
{
	int n_population, n_gene, n_constraint;
	char *out_solution, *out_population;
	DecodingFun decodingfun;

	n_population = Param->n_population;
	n_gene = Param->n_gene;
	n_constraint = Param->n_constraint;
	out_solution = Param->out_solution;
	out_population = Param->out_population;
	decodingfun = Param->decodingfun;

	writePopulation(out_population,Population,n_population,n_gene,n_constraint,decodingfun);
	writeSolution(t0,out_solution,best,n_gene,n_constraint,decodingfun);
	printElapsedTime(t0);
	free(t0);
	freePopulation(Population, n_population, n_gene, n_constraint);
	freeChrom(best,n_gene,n_constraint);
	free(Param);
	
	return;
}


void InitPopulation(RCGAParam *Param, Chromosome **Population)
{
	int i, j;

	for( i=1; i<=Param->n_population; i++ ){
		for( j=1; j<=Param->n_gene; j++ ) Population[i]->gene[j] = Rand();
		setFitness(Population[i],Param->n_gene,Param->n_constraint,Param->fitnessfun,Param->decodingfun);
	}
	sort(Population,Param->n_population,Param->n_gene,Param->n_constraint,Param->Pf);
}
