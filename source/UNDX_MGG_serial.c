#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "UNDX_MGG_serial.h"


void RCGA(RCGAParam *Param, Chromosome **Population, Chromosome *best, time_t *t0)
{
	int i, index;
	int n_population, n_gene, n_constraint, n_generation, output_intvl;
	double allowable_error;
	long t_limit;
	char *out_transition;
	DecodingFun decodingfun;
	int flg_printed = 0;

	n_population = Param->n_population;
	n_gene = Param->n_gene;
	n_constraint = Param->n_constraint;
	n_generation = Param->n_generation;
	output_intvl = Param->output_intvl;
	allowable_error = Param->allowable_error;
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
	if( (best->phi == 0 && best->f <= allowable_error) || getElapsedTime(t0) >= t_limit ) goto end;

	while(i < n_generation){
		i++;
		flg_printed = 0;
		MGG(Param,Population);
		index = findBest(Population,n_population);
		if( Population[index]->phi < best->phi 
			|| (Population[index]->phi == best->phi && Population[index]->f < best->f) )
			copyChrom(best,Population[index],n_gene,n_constraint);
		if(0 < output_intvl && i % output_intvl == 0){
			printTransition(t0,i,best->f,best->phi);
			writeTransition(t0,i,out_transition,best,n_gene,n_constraint,decodingfun);
			flg_printed = 1;
		}
		if( (best->phi == 0 && best->f <= allowable_error) || getElapsedTime(t0) >= t_limit ) goto end;
	}

end:
	if(0 < output_intvl && flg_printed == 0){
		printTransition(t0,i,best->f,best->phi);
		writeTransition(t0,i,out_transition,best,n_gene,n_constraint,decodingfun);
	}

	return;
}


void MGG(RCGAParam *Param, Chromosome **Population)
{
	/*
	Minimal Generation Gap Selection for UNDX
	See Hiroaki Kitano, "Genetic Algorithms 4", Sangyo-tosho, pp234, 2000 
	*/

	const int maxitr = 10;
	int *ip;
	Chromosome **f;
	int i, j;
	int n_children, n_gene, n_constraint, n_population;
	double Pf;

	n_children = Param->n_children;
	n_gene = Param->n_gene;
	n_constraint = Param->n_constraint;
	n_population = Param->n_population;
	Pf = Param->Pf;

	ip = allocIVector(3);
	f = allocPopulation(n_children+2,n_gene,n_constraint);

	// Picking up two parents from main population randomly
	for(i=1; i<=maxitr; i++){
		for(j=1; j<=2; j++) ip[j] = (int)( Rand() * (double)n_population ) + 1;
		if(ip[1] != ip[2]) break;
	}

	// Generating children
	for(i=1; i<=n_children; i++){
		// The third parent are randomly chosen each time.
		for(j=1; j<=maxitr; j++){
			ip[3] = (int)( Rand() * (double)n_population ) + 1;
			if(ip[3] != ip[1] && ip[3] != ip[2]) break;
		}
		getNewChild(Param,Population[ip[1]],Population[ip[2]],Population[ip[3]],f[i]);
	}
	// Adding main two parents into group of children 
	for(i=1; i<=2; i++) copyChrom(f[n_children+i],Population[ip[i]],n_gene,n_constraint);
	// Replacing main two parents with best individual and randomly chosen individual 
	sort(f,n_children+2,n_gene,n_constraint,Pf);
	copyChrom( Population[ip[1]], f[1], n_gene, n_constraint);
	copyChrom( Population[ip[2]], f[ (int)( Rand() * (double)(n_children+1) ) + 2 ], n_gene, n_constraint );
	sort(Population,n_population,n_gene,n_constraint,Pf);
	freePopulation(f,n_children+2,n_gene,n_constraint);
	freeIVector(ip,3);
	
	return;
}


void getNewChild(RCGAParam *Param, Chromosome *p1, Chromosome *p2, Chromosome *p3, Chromosome *c)
{
	/*
	Hander of UNDX function
	Checking wheather children are in search region 
	*/

	const int maxitr = 100;
	int flg;
	int i, j;

	for(i=1; i<=maxitr; i++){
		flg = 0;
		UNDX(Param,p1,p2,p3,c);
		for(j=1; j<=Param->n_gene ;j++)
			if(c->gene[j] < 0.0 || 1.0 < c->gene[j])
				flg=1;
		if(flg == 0) break;
	}
	if(flg == 1){
		for(j=1;j<=Param->n_gene;j++){
			if(1.0 < c->gene[j]){
				c->gene[j]=1.0;
			}else if(c->gene[j] < 0.0){
				c->gene[j]=0.0;
			}
		}
	}

	setFitness(c,Param->n_gene,Param->n_constraint,Param->fitnessfun,Param->decodingfun);

	return;
}


void UNDX(RCGAParam *Param, Chromosome *p1, Chromosome *p2, Chromosome *p3, Chromosome *c)
{
	/*
	Function of "Unimordal Distribution Crossover"
	p1 and p2 are main parents, and p3 is sub parent.
	See Hiroaki Kitano, "Genetic Algorithms 4", Sangyo-tosho, pp261, 2000
	*/

	const double alpha = 0.5, beta = 0.35; // Parameters for UNDX
	double d1, d2;
	double *v12, *v13, *e1, *t;
	double temp1, temp2;
	int i;
	int n_gene;

	n_gene = Param->n_gene;

	v12 = allocDVector(n_gene);
	v13 = allocDVector(n_gene);
	e1 = allocDVector(n_gene);
	t = allocDVector(n_gene);

	// v12 is a vector from p1 to p2, v13 is a vector from p1 to p3
	for(i=1; i<=n_gene; i++){
		v12[i] = p2->gene[i] - p1->gene[i];
		v13[i] = p3->gene[i] - p1->gene[i];
	}
	// temp1 is an inner product (v12,v13), temp2 is (v12,v12)
	temp1 = 0.0;
	temp2 = 0.0;
	for(i=1; i<=n_gene; i++){
		temp1 += v12[i] * v13[i];
		temp2 += v12[i] * v12[i];
	}
	d1 = sqrt(temp2); // d1 is distance between p1 and p2
	d2 = 0;
	if(temp2 == 0.0 || d1 == 0.0){
		for(i=1; i<=n_gene; i++){
			d2 += v13[i] * v13[i];
			e1[i] = 0.0;
		}
	}else{
		for(i=1; i<=n_gene; i++){
			d2 += ( temp1 / temp2 * v12[i] - v13[i] ) * ( temp1 / temp2 * v12[i] - v13[i] );
			e1[i] = v12[i] / d1;
		}
	}
	d2 = sqrt(d2); //d2 is distance between primary search line and p13

	for(i=1; i<=n_gene; i++)
		t[i] = NormalRand() * beta * d2 / sqrt( (double)n_gene );
	temp1 = 0.0; // temp1 is the inner product (t,e1)
	for(i=1; i<=n_gene; i++) temp1 += t[i] * e1[i];
	temp2 = alpha * d1 * NormalRand(); // temp2 is standard deviation for direction of e1
	for(i=1; i<=n_gene; i++)
		c->gene[i] = 0.5 * ( p1->gene[i] + p2->gene[i] ) + t[i] + ( temp2 - temp1 ) * e1[i];

	freeDVector(v12,n_gene);
	freeDVector(v13,n_gene);
	freeDVector(e1,n_gene);
	freeDVector(t,n_gene);

	return;
}


void RCGAInitial(int *argc, char ***argv, int seed, int n_gene, int n_generation, int n_population,
	int n_children, double allowable_error, long t_limit,
	int n_constraint, double Pf, int output_intvl, char *out_transition, char *out_solution, char *out_population,
	FitnessFun fitnessfun, DecodingFun decodingfun,
	RCGAParam **Param, Chromosome ***Population, Chromosome **best, time_t **t0)
{
	srand(seed);

	(*Param) = (RCGAParam *)malloc((size_t)sizeof(RCGAParam));
	(*Param)->n_gene = n_gene;
	(*Param)->n_generation = n_generation;
	(*Param)->n_population = n_population;
	(*Param)->n_children = n_children;
	(*Param)->allowable_error = allowable_error;
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
