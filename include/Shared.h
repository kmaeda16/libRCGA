#ifndef _Shared_h
#define _Shared_h

typedef void (*FitnessFun)(double *, int, int, double *, double *);
typedef void (*DecodingFun)(double *, int, double *);

typedef struct {
	double *gene;
	double f;
	double *g;
	double phi;
} Chromosome;

void setFitness(Chromosome *chrom, int n_gene, int n_constraint, FitnessFun fitnessfun, DecodingFun decodingfun);
void printTransition(time_t *t0, int generation, double f, double phi);
void writeTransition(time_t *t0, int generation, char *out_transition,
	Chromosome *chrom, int n_gene, int n_constraint, DecodingFun decodingfun);
void writePopulation(char *out_population, Chromosome **Population,
	int n_population, int n_gene, int n_constraint, DecodingFun decodingfun);
void writeSolution(time_t *t0, char *out_solution,
	Chromosome *chrom, int n_gene, int n_constraint, DecodingFun decodingfun);
void sort(Chromosome **Population, int n_population, int n_gene, int n_constraint, double Pf);
void HeapSort_PfZero(Chromosome **Population, int n_population, int n_gene, int n_constraint);
void HeapSort_PfOne(Chromosome **Population, int n_population, int n_gene, int n_constraint);
void StochasticRankingSort(Chromosome **Population, int n_population, int n_gene, int n_constraint, double Pf);
int findBest(Chromosome **Population, int n_population);
double Rand();
double NormalRand();
void printElapsedTime(time_t *t0);
long getElapsedTime(time_t *t0);
Chromosome **allocPopulation(int n_population, int n_gene, int n_constraint);
void freePopulation(Chromosome **Population, int n_population, int n_gene, int n_constraint);
Chromosome *allocChrom(int n_gene, int n_constraint);
void freeChrom(Chromosome *chrom, int n_gene, int n_constraint);
void copyChrom(Chromosome *to, Chromosome *from, int n_gene, int n_constraint);
double *allocDVector(int n);
void freeDVector(double *v, int n);
int *allocIVector(int n);
void freeIVector(int *v, int n);
void dswap(double *a, double *b);
void iswap(int *a, int *b);
void handleError(char *errormsg);

#endif
