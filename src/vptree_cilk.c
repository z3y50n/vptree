#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <time.h>
#include <sys/time.h>

#include <cilk/cilk.h>

#include "../inc/vptree.h"
#include "../inc/quickselect.h"

vptree * recbuildvp(double *X, int *index, int n, int d);
double * findDistance(double *X, int n, int d);
double median(double *distance, int n);


vptree * buildvp(double *X, int n, int d)
{
	vptree *T = malloc(sizeof(vptree));
	// Create array of indexes
	int *index = malloc(n*sizeof(int));
	for(int i=0;i<n;i++)
	{
		*(index+i) = i;
	}

	T = recbuildvp(X, index, n, d);
	
	free(index);
	return T;
}

vptree * recbuildvp(double *X, int *index, int n, int d)
{
	vptree *T = NULL;
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
		return T;	
	}

	T = calloc(1,sizeof(vptree));	
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
		return T;
	}

	
	//Calculate Distances
	staticdistance = malloc((n-1)*sizeof(double));
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
	innerindex = malloc(innersize*sizeof(int));
	inner = malloc(innersize*d*sizeof(double));
	outerindex = malloc(outersize*sizeof(int));
	outer = malloc(outersize*d*sizeof(double));
	//Create Inner and Outer arrays
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


	//Print Nodes of VP Tree
	// printf("\t\t=======%d=======\n",n );
	// printf("\t\tidx: %d\n",T->idx );
	// printf("\t\tmedian: %f\n", T->md );
	// printf("\t\tvp: %f %f \n",T->vp[0],T->vp[1]);

	T->inner = cilk_spawn recbuildvp(inner,innerindex,innersize,d);
	T->outer = cilk_spawn recbuildvp(outer,outerindex,outersize,d);
	cilk_sync;

	free(distance);
	free(staticdistance);
	free(inner);
	free(outer);
	free(innerindex);
	free(outerindex);

	return T;
}


double * findDistance(double *X, int n, int d)
{
	double *distance = calloc((n-1),sizeof(double));
	if(n>1000)
	{
		cilk_for (int i=0;i<n-1;i++)
		{
			for(int j=0;j<d;j++)
			{
				*(distance+i) += (*(X+n*j+i) - *(X+n*j+n-1))*(*(X+n*j+i) - *(X+n*j+n-1));
			}
			*(distance+i) = sqrt(*(distance+i));
		}
	}
	else
	{
		for (int i=0;i<n-1;i++)
		{
			for(int j=0;j<d;j++)
			{
				*(distance+i) += (*(X+n*j+i) - *(X+n*j+n-1))*(*(X+n*j+i) - *(X+n*j+n-1));
			}
			*(distance+i) = sqrt(*(distance+i));
		}
	}
	
	return distance;
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