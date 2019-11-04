#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <time.h>
#include <sys/time.h>

#include <pthread.h>

#include "../inc/vptree.h"
#include "../inc/quickselect.h"

#define NUM_THREADS 10
#define THRESHOLD 1000
#define MAX_THREAD 20

int activethreads = 0;

void * recbuildvp(void *arg);
double * findDistance(double *X, int n, int d);
void *threadDistance(void * args);
double median(double *distance, int n);


vptree * buildvp(double *X, int n, int d)
{
	void *T = malloc(sizeof(vptree));
	build *arg = calloc(1,sizeof(build));

	// Create array of indexes
	int *index = malloc(n*sizeof(int));
	for(int i=0;i<n;i++)
	{
		*(index+i) = i;
	}
	arg->X = malloc(n*d*sizeof(double));
	arg->index = malloc(n*sizeof(int));
	arg->X = X;
	arg->index = index;
	arg->n = n;
	arg->d = d;

	T = recbuildvp((void*)arg);
	
	free(index);
	return (vptree *)T;
}

void * recbuildvp(void *arg)
{
	vptree *T = NULL;
	build *args = (build *) arg;
	void *res;
	int d = args->d;
	int n = args->n;
	int *index = args->index;
	double *X = args->X;

	build *arginner = calloc(1,sizeof(build));
	build *argouter = calloc(1,sizeof(build));

	double *distance;
	double *staticdistance;
	int *innerindex;
	int *outerindex;
	int innersize = 0;
	int outersize = 0;
	double *inner;
	double *outer;	

	if(n==0)
	{
		return (void *)T;	
	}

	T = calloc(1,sizeof(vptree));	
	res = calloc(1,sizeof(vptree));
	//Find Vantage Point
	T->vp = malloc(d*sizeof(double));
	for(int i=0;i<d;i++)
	{
		T->vp[i]=*(X+n*i+n-1);
	}

	//Index of Vantage Point
	T->idx = *(index+n-1);
	
	if(n==1)
	{
		return (void *)T;
	}

	//Calculate Distances
	staticdistance = calloc(n-1,sizeof(double));
	distance = findDistance(X,n,d);
	for(int i=0;i<n-1;i++)
	{
		staticdistance[i] = distance[i];
	}


	//Find Median
	T->md = median(distance,n);

	//Calculate size of inner/outer
	for(int i=0;i<n-1;i++)
	{
		if(staticdistance[i]<=T->md)
		{
			innersize++;
		}
		else
		{
			outersize++;
		}
	}
	
	//Create Inner and Outer arrays
	innerindex = malloc(innersize*sizeof(int));
	inner = malloc(innersize*d*sizeof(double));
	outerindex = malloc(outersize*sizeof(int));
	outer = malloc(outersize*d*sizeof(double));
	
	int in = 0;
	int out = 0;
	for(int i=0;i<n-1;i++)
	{
		if(staticdistance[i]<=T->md)
		{

			*(innerindex+in) = *(index+i);
			for(int j=0;j<d;j++)
			{
				*(inner+in+innersize*j) = *(X+i+n*j);
			}
			in++;
		}
		else
		{
			*(outerindex+out) = *(index+i);
			for(int j=0;j<d;j++)
			{
				*(outer+out+outersize*j) = *(X+i+n*j);
			}
			out++;
		}
	}


	arginner->n = innersize;
	arginner->d = d;
	arginner->index = innerindex;
	arginner->X = inner;

	argouter->n = outersize;
	argouter->d = d;
	argouter->index = outerindex;
	argouter->X = outer;

	//Print Nodes of VP Tree
	// printf("\t\t=======%d=======\n",n );
	// printf("\t\tidx: %d\n",T->idx );
	// printf("\t\tmedian: %f\n", T->md );
	// printf("\t\tvp: %f %f \n",T->vp[0],T->vp[1]);

	pthread_t t_in;
	if(activethreads>MAX_THREAD || n<10)
	{
		T->inner = (vptree *)recbuildvp((void *)arginner);
		T->outer = (vptree *)recbuildvp((void *)argouter);
		free(distance);
		free(staticdistance);
		free(inner);
		free(outer);
		free(innerindex);
		free(outerindex);
		return (void *)T;
	}
	activethreads++;
	pthread_create(&t_in,NULL,recbuildvp,(void *)arginner);
	T->outer = (vptree *)recbuildvp((void *)argouter);
	pthread_join(t_in,&res);
	activethreads--;
	T->inner = (vptree *) res;
	
	free(distance);
	free(staticdistance);
	free(inner);
	free(outer);
	free(innerindex);
	free(outerindex);

	return (void *)T;
}

double * findDistance(double *X, int n, int d)
{
	double *distance = calloc((n-1),sizeof(double));
	if(n<THRESHOLD)
	{
		for(int i=0;i<n-1;i++)
		{
			for(int j=0;j<d;j++)
			{
				*(distance+i) += (*(X+n*j+i) - *(X+n*j+n-1))*(*(X+n*j+i) - *(X+n*j+n-1));
			}
			*(distance+i) = sqrt(*(distance+i));
		}
		return distance;
	}

	ForDistance *arg;
	arg = calloc(NUM_THREADS,sizeof(ForDistance));
	pthread_t threads[NUM_THREADS];
	
	//Initialize arg
	for(int i=0;i<NUM_THREADS;i++)
	{
		arg[i].n = n;
		arg[i].d = d;
		arg[i].start = ((n-1)/NUM_THREADS)*i;
		arg[i].stop = arg[i].start + (n-1)/NUM_THREADS ;
		if(i==NUM_THREADS-1)
		{
			arg[i].stop = n-1;
		}
		arg[i].distance = distance;
		arg[i].X = X;
		pthread_create(&threads[i],NULL,threadDistance, &arg[i]);
	}

	for(int i=0;i<NUM_THREADS;i++)
	{
		pthread_join(threads[i],NULL);
	}
	return arg[0].distance;
	
}

void *threadDistance(void * args)
{
	ForDistance *arg = (ForDistance *) args;
	for(int i=arg->start;i<arg->stop;i++)
	{
		for(int j=0;j<arg->d;j++)
		{
			arg->distance[i] += (arg->X[arg->n*j+i] - arg->X[arg->n*j+arg->n-1])*(arg->X[arg->n*j+i] - arg->X[arg->n*j+arg->n-1]);
		}
		arg->distance[i] = sqrt(arg->distance[i]);	
	}
}
double median(double *distance, int n)
{
	double md;
	if((n-1)%2!=0)
	{
		md = ksmallest(distance,n-1,(n-1)/2+1);
	}
	else
	{
		md = (ksmallest(distance,n-1,(n-1)/2+1)+ksmallest(distance,n-1,(n-1)/2))/2;
	}
	return md;
}

vptree * getInner(vptree * T)
{
	return T->inner;
}
vptree * getOuter(vptree * T)
{
	return T->outer;
}
double getMD(vptree * T)
{
	return T->md;
}
double * getVP(vptree * T)
{
	return T->vp;
}
int getIDX(vptree * T)
{
	return T->idx;
}