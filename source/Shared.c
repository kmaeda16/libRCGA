#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <string.h>
#include "Shared.h"


void setFitness(Chromosome *chrom, int n_gene, int n_constraint, FitnessFun fitnessfun, DecodingFun decodingfun)
{
	double *x;
	int i;

	x = allocDVector(n_gene);

	decodingfun(chrom->gene,n_gene,x);
	fitnessfun(x,n_gene,n_constraint,&chrom->f,chrom->g);

	chrom->phi = 0;
	for(i=1; i<=n_constraint; i++) chrom->phi += fmax(0,chrom->g[i]) * fmax(0,chrom->g[i]);

	freeDVector(x,n_gene);

	return;
}


void printTransition(time_t *t0, int generation, double f, double phi)
{
	long int t;

	t = getElapsedTime(t0);
	printf("Elapsed Time = %ld, Generation = %d, f = %e, phi = %e\n",t,generation,f,phi);

	return;
}


void writeTransition(time_t *t0, int generation, char *out_transition,
	Chromosome *chrom, int n_gene, int n_constraint, DecodingFun decodingfun)
{
	double *x;
	FILE *out;
	int i;
	long int t;

	if( strcmp("None",out_transition) == 0 || strcmp("none",out_transition) == 0)
		return;

	if(generation == 1){
		if( ( out = fopen(out_transition,"w") ) == NULL ){
			fprintf(stderr,"cannot open %s!\n",out_transition);
			return;
		}
		fprintf(out,"Time\tGeneration\tf\tphi\t");
		for(i=1; i<=n_gene      ; i++) fprintf(out,"x[%d]\t",i);
		for(i=1; i<=n_constraint; i++) fprintf(out,"g[%d]\t",i);
		fprintf(out,"\n");
	}else{
		if( ( out = fopen(out_transition,"a") ) == NULL ){
			fprintf(stderr,"cannot open %s!\n",out_transition);
			return;
		}
	}

	t = getElapsedTime(t0);
	x = allocDVector(n_gene);
	decodingfun(chrom->gene,n_gene,x);

	fprintf(out,"%ld\t%d\t%e\t%e\t",t,generation,chrom->f,chrom->phi);
	for(i=1; i<=n_gene       ; i++) fprintf(out,"%e\t",x[i]);
	for(i=1; i<=n_constraint ; i++) fprintf(out,"%e\t",chrom->g[i]);
	fprintf(out,"\n");

	freeDVector(x,n_gene);
	fclose(out);

	return;
}


void writePopulation(char *out_population, Chromosome **Population,
	int n_population, int n_gene, int n_constraint, DecodingFun decodingfun)
{
	double *x;
	FILE *out;
	int i, j;

	if( strcmp("None",out_population) == 0 || strcmp("none",out_population) == 0)
		return;

	if( ( out = fopen(out_population,"w") ) == NULL ){
		fprintf(stderr,"cannot open %s!\n",out_population);
		return;
	}

	x = allocDVector(n_gene);
	
	fprintf(out,"No\tf\tphi\t");
	for( j=1; j<=n_gene      ; j++ ) fprintf(out,"x[%d]\t",j);
	for( j=1; j<=n_constraint; j++ ) fprintf(out,"g[%d]\t",j);
	fprintf(out,"\n");
	
	for( i=1; i<=n_population; i++ ){
		decodingfun(Population[i]->gene,n_gene,x);
		fprintf(out,"%d\t%e\t%e\t",i,Population[i]->f,Population[i]->phi);
		for( j=1; j<=n_gene      ; j++ ) fprintf(out,"%e\t",x[j]);
		for( j=1; j<=n_constraint; j++ ) fprintf(out,"%e\t",Population[i]->g[j]);
		fprintf(out,"\n");
	}

	freeDVector(x,n_gene);
	fclose(out);

	return;
}


void writeSolution(time_t *t0, char *out_solution,
	Chromosome *chrom, int n_gene, int n_constraint, DecodingFun decodingfun)
{
	double *x;
	FILE *out;
	int i;
	long int t;

	if( strcmp("None",out_solution) == 0 || strcmp("none",out_solution) == 0)
		return;

	if( ( out = fopen(out_solution,"w") ) == NULL ){
		fprintf(stderr,"cannot open %s!\n",out_solution);
		return;
	}

	t = getElapsedTime(t0);
	x = allocDVector(n_gene);
	decodingfun(chrom->gene,n_gene,x);

	fprintf(out,"Time\tf\tphi\t");
	for( i=1; i<=n_gene      ; i++ ) fprintf(out,"x[%d]\t",i);
	for( i=1; i<=n_constraint; i++ ) fprintf(out,"g[%d]\t",i);
	fprintf(out,"\n");

	fprintf(out,"%ld\t%e\t%e\t",t,chrom->f,chrom->phi);
	for( i=1; i<=n_gene      ; i++ ) fprintf(out,"%e\t",x[i]);
	for( i=1; i<=n_constraint; i++ ) fprintf(out,"%e\t",chrom->g[i]);
	fprintf(out,"\n");

	freeDVector(x,n_gene);
	fclose(out);

	return;
}


void sort(Chromosome **Population, int n_population, int n_gene, int n_constraint, double Pf)
{
	/*
	If Pf <= 0 or 1 <= Pf, the stochastic ranking does not show any stochasticity.
	In this case, HeapSort_PfZero or HeapSort_PfOne are recommended because they
	are faster and more memory efficient than StochasticRankingSort.
	If n_constraint <= 0, only f is used for ranking.
	*/
	
	if(1 <= Pf || n_constraint <= 0){ 
		// f always dominates phi (only f is evaluated)
		HeapSort_PfOne(Population,n_population,n_gene,n_constraint);
	} else if(Pf <= 0){
		// phi always dominates f (firstly, phi is evaluated, then f is evaluated)
		HeapSort_PfZero(Population,n_population,n_gene,n_constraint);
	}else{
		StochasticRankingSort(Population,n_population,n_gene,n_constraint,Pf);
	}

}


void HeapSort_PfZero(Chromosome **Population, int n_population, int n_gene, int n_constraint)
{
	int l,j,ir,i;
	Chromosome *tmp;

	tmp = allocChrom(n_gene,n_constraint);

	l=(n_population >> 1)+1;
	ir=n_population;
	for (;;) {
		if (l > 1)
			copyChrom(tmp,Population[--l],n_gene,n_constraint);
		else {
			copyChrom(tmp,Population[ir],n_gene,n_constraint);
			copyChrom(Population[ir],Population[1],n_gene,n_constraint);
			if (--ir == 1) {
				copyChrom(Population[1],tmp,n_gene,n_constraint);
				return;
			}
		}
		i=l;
		j=l << 1;
		while (j <= ir) {
			if (j < ir && (Population[j]->phi < Population[j+1]->phi || (Population[j]->phi == Population[j+1]->phi && Population[j]->f < Population[j+1]->f)) ) ++j;
			if ( tmp->phi < Population[j]->phi || (tmp->phi == Population[j]->phi && tmp->f < Population[j]->f) ) {
				copyChrom(Population[i],Population[j],n_gene,n_constraint);
				j += (i=j);
			}
			else j=ir+1;
		}
		copyChrom(Population[i],tmp,n_gene,n_constraint);
	}

	freeChrom(tmp,n_gene,n_constraint);
}


void HeapSort_PfOne(Chromosome **Population, int n_population, int n_gene, int n_constraint)
{
	int l,j,ir,i;
	Chromosome *tmp;

	tmp = allocChrom(n_gene,n_constraint);

	l=(n_population >> 1)+1;
	ir=n_population;
	for (;;) {
		if (l > 1)
			copyChrom(tmp,Population[--l],n_gene,n_constraint);
		else {
			copyChrom(tmp,Population[ir],n_gene,n_constraint);
			copyChrom(Population[ir],Population[1],n_gene,n_constraint);
			if (--ir == 1) {
				copyChrom(Population[1],tmp,n_gene,n_constraint);
				return;
			}
		}
		i=l;
		j=l << 1;
		while (j <= ir) {
			if (j < ir && Population[j]->f < Population[j+1]->f ) ++j;
			if ( tmp->f < Population[j]->f ) {
				copyChrom(Population[i],Population[j],n_gene,n_constraint);
				j += (i=j);
			}
			else j=ir+1;
		}
		copyChrom(Population[i],tmp,n_gene,n_constraint);
	}

	freeChrom(tmp,n_gene,n_constraint);
}


void StochasticRankingSort(Chromosome **Population, int n_population, int n_gene, int n_constraint, double Pf)
{
	Chromosome **pop_tmp;
	int *index_tmp;
	double *f_tmp, *phi_tmp;
	double u;
	int i, j;
	int index1, index2;
	int flg;

	index_tmp = allocIVector(n_population);
	f_tmp = allocDVector(n_population);
	phi_tmp = allocDVector(n_population);

	for(i=1; i<=n_population; i++){
		index_tmp[i] = i;
		f_tmp[i] = Population[i]->f;
		phi_tmp[i] = Population[i]->phi;
	}
	
	for(i=1; i<=n_population-1; i++){
		index1 = (int)(Rand()*(double)n_population) + 1;
		index2 = (int)(Rand()*(double)n_population) + 1;
		iswap(&index_tmp[index1],&index_tmp[index2]);
		dswap(&f_tmp[index1],&f_tmp[index2]);
		dswap(&phi_tmp[index1],&phi_tmp[index2]);
	}
	
	for(i=1; i<=n_population; i++){
		flg = 0;
		for(j=1; j<=n_population-1; j++){
			u = Rand();
			if( phi_tmp[j] == phi_tmp[j+1] || u < Pf ){
				if( f_tmp[j] > f_tmp[j+1] ){
					iswap(&index_tmp[j],&index_tmp[j+1]);
					dswap(&f_tmp[j],&f_tmp[j+1]);
					dswap(&phi_tmp[j],&phi_tmp[j+1]);
					flg = 1;
				}
			}else{
				if( phi_tmp[j] >  phi_tmp[j+1] ){
					iswap(&index_tmp[j],&index_tmp[j+1]);
					dswap(&f_tmp[j],&f_tmp[j+1]);
					dswap(&phi_tmp[j],&phi_tmp[j+1]);
					flg = 1;
				}
			}
		}
		if(flg == 0) break;
	}
	
	pop_tmp = allocPopulation(n_population,n_gene,n_constraint);

	for(i=1; i<=n_population; i++) copyChrom(pop_tmp[i],Population[i],n_gene,n_constraint);
	for(i=1; i<=n_population; i++) copyChrom(Population[i],pop_tmp[index_tmp[i]],n_gene,n_constraint);

	freePopulation(pop_tmp,n_population,n_gene,n_constraint);
	freeIVector(index_tmp,n_population);
	freeDVector(f_tmp,n_population);
	freeDVector(phi_tmp,n_population);

}


int findBest(Chromosome **Population, int n_population)
{
	double f = DBL_MAX, phi = DBL_MAX;
	int i, index = 0;

	for(i=n_population; 1<=i; i--){
		if( Population[i]->phi < phi
			|| (Population[i]->phi == phi && Population[i]->f < f) ){
			f = Population[i]->f;
			phi = Population[i]->phi;
			index = i;
		}
	}

	return(index);

}


double Rand()
{
	return ( (double)rand() + 1.0 ) / ( (double)RAND_MAX + 2.0 );
}


double NormalRand()
{
	static int iset = 0;
	static double gset;
	double fac, rsq, v1, v2;

	if(iset == 0){
		do{
			v1 = 2.0 * Rand() - 1.0;
			v2 = 2.0 * Rand() - 1.0;
			rsq = v1 * v1 + v2 * v2;
		}while(rsq >= 1.0 || rsq == 0.0);
		fac = sqrt( -2.0 * log(rsq) / rsq );
		gset = v1 * fac;
		iset = 1;
		return v2 * fac;
	}else{
		iset = 0;
		return gset;
	}
}


void printElapsedTime(time_t *t0)
{
	long t, h, m, s;

	t = getElapsedTime(t0);

	h = t / 3600;
	m = ( ( t % 3600 ) / 60 );
	s = ( ( t % 3600 ) % 60 );

	printf("Elapsed Time: %ldhour %ldmin %ldsec\n",h,m,s);
}


long getElapsedTime(time_t *t0)
{
	return( (long)(time(NULL) - (*t0)) );
}


Chromosome **allocPopulation(int n_population, int n_gene, int n_constraint)
{
	Chromosome **Population;
	int i;

	Population = (Chromosome **) malloc( (size_t)( n_population * sizeof(Chromosome *) ) );
	for(i=1; i<=n_population; i++) Population[i-1] = allocChrom(n_gene,n_constraint);

	return(Population-1);
}


void freePopulation(Chromosome **Population, int n_population, int n_gene, int n_constraint)
{
	int i;

	for(i=1; i<=n_population; i++) freeChrom(Population[i],n_gene,n_constraint);
	free(Population+1);

	return;
}


Chromosome *allocChrom(int n_gene, int n_constraint)
{
	Chromosome *chrom;

	chrom = (Chromosome *) malloc( (size_t) sizeof(Chromosome) );
	chrom->gene = allocDVector(n_gene);
	chrom->g = allocDVector(n_constraint);

	return(chrom);
}


void freeChrom(Chromosome *chrom, int n_gene, int n_constraint)
{
	freeDVector(chrom->gene,n_gene);
	freeDVector(chrom->g,n_constraint);
	free(chrom);

	return;
}


void copyChrom(Chromosome *to, Chromosome *from, int n_gene, int n_constraint)
{
	int i;

	for(i=1; i<=n_gene; i++) to->gene[i] = from->gene[i];
	to->f = from->f;
	for(i=1; i<=n_constraint; i++) to->g[i] = from->g[i];
	to->phi = from->phi;

	return;
}


double *allocDVector(int n)
{
	double *v;

	v = (double *) malloc( (size_t) ( n * sizeof(double) ) );
	if(!v) handleError("Error in allocDVector.");

	return(v-1);
}


void freeDVector(double *v, int n)
{
	free(v+1);
}


int *allocIVector(int n)
{
	int *v;

	v = (int *) malloc( (size_t) ( n * sizeof(int) ) );
	if(!v) handleError("Error in allocIVector.");

	return(v-1);
}


void freeIVector(int *v, int n)
{
	free(v+1);
}


void dswap(double *a, double *b)
{
	double c;

	c = (*a);
	(*a) = (*b);
	(*b) = c;

}


void iswap(int *a, int *b)
{
	int c;

	c = (*a);
	(*a) = (*b);
	(*b) = c;

}


void handleError(char *errormsg)
{
	fprintf(stderr,"%s\n",errormsg);
	fprintf(stderr,"Exit!\n");
	exit(1);
}
