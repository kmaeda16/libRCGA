#define _USE_MATH_DEFINES
#define INIT_VALUE (-1e+10)
#define DISMISSING_SIGNAL (-1e+20)
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include "UNDX_MGG_parallel.h"


void RCGA(RCGAParam *Param, Chromosome **Population, Chromosome *best, time_t *t0)
{
	int i, index;
	int n_population, n_gene, n_constraint, n_generation, output_intvl;
	double allowable_error;
	long t_limit;
	char *out_transition;
	DecodingFun decodingfun;
	int flg_printed = 0;

	if(Param->my_rank > 0){
		calcFitness4Master(Param);
		return;
	}

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
	dismissWorkers(Param);

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
	requestFitnessCalc(Param,f,n_children);
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

	//setFitness(c,Param->n_gene,Param->n_constraint,Param->fitnessfun,Param->decodingfun);

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
	int my_rank, n_ps;

	MPI_Init(argc, argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &n_ps);

	if (n_ps < 2) handleError("Error in RCGAInitial: # processes must be > 1!");
	
	srand(seed);

	(*Param) = (RCGAParam *)malloc((size_t)sizeof(RCGAParam));
	(*Param)->my_rank = my_rank;
	(*Param)->n_ps = n_ps;
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
	
	if (my_rank == 0) {
		(*Population) = allocPopulation(n_population, n_gene, n_constraint);
		(*best) = allocChrom(n_gene,n_constraint);
		(*t0) = (time_t *)malloc((size_t)sizeof(time_t));
		(**t0) = time(NULL);
	}

	return;
}


void RCGADeInitial(RCGAParam *Param, Chromosome **Population, Chromosome *best, time_t *t0)
{
	int my_rank;
	int n_population, n_gene, n_constraint;
	char *out_solution, *out_population;
	DecodingFun decodingfun;

	my_rank = Param->my_rank;
	n_population = Param->n_population;
	n_gene = Param->n_gene;
	n_constraint = Param->n_constraint;
	out_solution = Param->out_solution;
	out_population = Param->out_population;
	decodingfun = Param->decodingfun;

	if(my_rank == 0){
		writePopulation(out_population,Population,n_population,n_gene,n_constraint,decodingfun);
		writeSolution(t0,out_solution,best,n_gene,n_constraint,decodingfun);
		printElapsedTime(t0);
		free(t0);
		freePopulation(Population, n_population, n_gene, n_constraint);
		freeChrom(best,n_gene,n_constraint);
	}
	free(Param);

	MPI_Finalize();
	return;
}


void InitPopulation(RCGAParam *Param, Chromosome **Population)
{
	int i, j;

	for( i=1; i<=Param->n_population; i++ ){
		for( j=1; j<=Param->n_gene; j++ ) Population[i]->gene[j] = Rand();
		// setFitness(Population[i],Param->n_gene,Param->n_constraint,Param->fitnessfun,Param->decodingfun);
	}
	requestFitnessCalc(Param,Population, Param->n_population);
	sort(Population,Param->n_population,Param->n_gene,Param->n_constraint,Param->Pf);
}


void requestFitnessCalc(RCGAParam *Param, Chromosome **Group, int n_group)
{
	MPI_Status status;
	double *buffer_send, *buffer_recv;
	int n_ps, n_gene, n_constraint, n_population, n_children, n_assignment;
	int *index;
	int count;
	int i, j, k;

	if(n_group == 0) return;

	n_ps = Param->n_ps;
	n_gene = Param->n_gene;
	n_constraint = Param->n_constraint;
	n_population = Param->n_population;
	n_children = Param->n_children;

	if(n_population > n_children){
		n_assignment = (int) ceil( (double)n_population / ((double)n_ps-1.0) );
	}else{
		n_assignment = (int) ceil( (double)n_children / ((double)n_ps-1.0) );
	}

#ifdef DEBUG /*-------------------------------------------------------------------------------------------------------*/
	for(j=1; j<=n_group; j++){
		printf("Group[%d]: %e\t%e:\t",j,Group[j]->f,Group[j]->phi);
		for(k=1; k<=n_gene; k++) printf("%e\t", Group[j]->gene[k]);
		for(k=1; k<=n_constraint; k++) printf("%e\t", Group[j]->g[k]);
		printf("\n");
	}
#endif /*-------------------------------------------------------------------------------------------------------------*/
	
	// index[j] has the index of process that will evaluates the j-th group
	index = allocIVector(n_group);
	for(j=1; j<=n_group; j++) index[j] = (j-1) % (n_ps-1) + 1;

	// Send some groups to workers
	buffer_send = allocDVector(n_assignment*n_gene);

	for(i=1; i<n_ps; i++){

		//	Subgroup -> buffer_send
		count = 0;
		for(j=1; j<=n_group; j++){
			if(index[j] == i){
				for(k=1; k<=n_gene; k++) buffer_send[++count] = Group[j]->gene[k];
			}
		}

#ifdef DEBUG /*-------------------------------------------------------------------------------------------------------*/
		for(j=1; j<=count; j++) printf("[%d] [0->%d]\tbuffer_send[%d]:\t%e\n",Param->my_rank,i,j,buffer_send[j]);
#endif /*-------------------------------------------------------------------------------------------------------------*/

		// Send
		MPI_Send(&buffer_send[1], count, MPI_DOUBLE, i, 0,MPI_COMM_WORLD);

	}
	freeDVector(buffer_send,n_assignment*n_gene);

	// Receive some groups from workers
	buffer_recv = allocDVector(n_assignment*(n_constraint+2));
	for(i=1; i<n_ps; i++){

		// Receive
		MPI_Recv(&buffer_recv[1], n_assignment*(n_constraint+2), MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);

#ifdef DEBUG /*-------------------------------------------------------------------------------------------------------*/
		for(j=1; j<=n_assignment*(n_constraint+2); j++) printf("[%d] [%d->0]\tbuffer_recv[%d]:\t%e\n",Param->my_rank,i,j,buffer_recv[j]);
#endif /*-------------------------------------------------------------------------------------------------------------*/

		// buffer_recv -> Group
		count = 0;
		for(j=1; j<=n_group; j++){
			if(index[j] == i){
				Group[j]->f = buffer_recv[++count];
				for(k=1; k<=n_constraint; k++) Group[j]->g[k] = buffer_recv[++count];
				Group[j]->phi = buffer_recv[++count];
			}
		}

	}
	freeDVector(buffer_recv,n_assignment*(n_constraint+2));

	freeIVector(index,n_group);

#ifdef DEBUG /*-------------------------------------------------------------------------------------------------------*/
	for(j=1; j<=n_group; j++){
		printf("Group[%d]: %e\t%e:\t",j,Group[j]->f,Group[j]->phi);
		for(k=1; k<=n_gene; k++) printf("%e\t", Group[j]->gene[k]);
		for(k=1; k<=n_constraint; k++) printf("%e\t", Group[j]->g[k]);
		printf("\n");
	}
#endif /*-------------------------------------------------------------------------------------------------------------*/

}


void calcFitness4Master(RCGAParam *Param)
{
	MPI_Status status;
	Chromosome *temp;
	int count, count1, count2;
	int j, k;
	double *buffer_recv, *buffer_send;
	int n_ps, n_population, n_children, n_gene, n_constraint, n_assignment, n_evalulated;

	n_ps = Param->n_ps;
	n_population = Param->n_population;
	n_children = Param->n_children;
	n_gene = Param->n_gene;
	n_constraint = Param->n_constraint;
	
	if(n_population > n_children){
		n_assignment = (int) ceil( (double)n_population / ((double)n_ps-1.0) );
	}else{
		n_assignment = (int) ceil( (double)n_children / ((double)n_ps-1.0) );
	}

	buffer_recv = allocDVector(n_assignment*n_gene);
	buffer_send = allocDVector(n_assignment*(n_constraint+2));
	temp = allocChrom(n_gene,n_constraint);

	while(1){
		
		// Initialization
		count = 0;
		for(j=1; j<=n_assignment; j++){
			for(k=1; k<=n_gene; k++) buffer_recv[++count] = INIT_VALUE;
		}
		
		// Receive
		MPI_Recv(&buffer_recv[1], n_assignment*n_gene, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
		
#ifdef DEBUG /*-------------------------------------------------------------------------------------------------------*/
		for(j=1; j<=n_assignment*n_gene; j++) printf("[%d] [0->%d]\tbuffer_recv[%d]:\t%e\n",Param->my_rank,Param->my_rank,j,buffer_recv[j]);
#endif /*-------------------------------------------------------------------------------------------------------------*/

		count = 0;
		for(j=1; j<=n_assignment; j++){
			for(k=1; k<=n_constraint+2; k++) buffer_send[++count] = INIT_VALUE;
		}

		count1 = 0;
		count2 = 0;
		n_evalulated = 0;
		for(j=1; j<=n_assignment; j++){

			for(k=1; k<=n_gene; k++) temp->gene[k] = buffer_recv[++count1];

			// If dismissing signal, break from while(1)
			if(temp->gene[1] == DISMISSING_SIGNAL) goto end;

			// Fitness calculation
			if(temp->gene[1] == INIT_VALUE) break;
			setFitness(temp,n_gene,n_constraint,Param->fitnessfun,Param->decodingfun);
			n_evalulated++;

			buffer_send[++count2] = temp->f;
			for(k=1; k<=n_constraint; k++) buffer_send[++count2] = temp->g[k];
			buffer_send[++count2] = temp->phi;
		}

		// Send
#ifdef DEBUG /*-------------------------------------------------------------------------------------------------------*/
		for(j=1; j<=count2; j++) printf("[%d] [%d->0]\tbuffer_send[%d]:\t%e\n",Param->my_rank,Param->my_rank,j,buffer_send[j]);
#endif /*-------------------------------------------------------------------------------------------------------------*/

		MPI_Send(&buffer_send[1], n_evalulated*(n_constraint+2), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		
	}

end:
	freeDVector(buffer_recv,n_assignment*n_gene);
	freeDVector(buffer_send,n_assignment*(n_constraint+2));
	freeChrom(temp,n_gene,n_constraint);

	return;
}


void dismissWorkers(RCGAParam *Param)
{
	double dismissing_signal;
	int i;

	dismissing_signal = DISMISSING_SIGNAL;
	for(i=1; i<Param->n_ps; i++)
		MPI_Send(&dismissing_signal, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);

	return;
}
