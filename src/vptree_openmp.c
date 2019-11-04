#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <omp.h>

#include <time.h>
#include <sys/time.h>

#include "../inc/vptree.h"
#include "../inc/quickselect.h"

#define N 100000
#define D 200

vptree * recbuildvp(double *X, int *index, int n, int d);
double * findDistance(double *X, int n, int d);
double median(double *distance, int n);


int main()
{

	// double X[20] = {0.840188, 0.394383, 0.783099, 0.798440, 0.911647, 0.197551, 0.335223, 0.768230, 0.277775, 0.553970,
 //      0.477397, 0.628871, 0.364784, 0.513401, 0.952230, 0.916195, 0.635712, 0.717297, 0.141603, 0.606969};
	srand(time(NULL));
	double *X = malloc(N * D * sizeof(double));
	for(int i=0;i<D;i++)
	{
		for(int j=0;j<N;j++)
		{
			*(X+i*N+j) = (double)(rand()%10+1);
			//printf("%f\n",*(X+i*N+j) );
		}
	}
	struct timeval start,end;
    double time=0;
    gettimeofday(&start,NULL);
	buildvp(X,N,D);
	gettimeofday(&end,NULL);
	time = (double)((end.tv_usec-start.tv_usec)/1.0e6 + end.tv_sec - start.tv_sec);
	printf("Time: %f\n",time );
	free(X);
	return 0;
}

vptree * buildvp(double *X, int n, int d)
{
	vptree *T = malloc(sizeof(vptree));
	// Create array of indexes
	int *index = malloc(n*sizeof(int));
	for(int i=0;i<n;i++)
	{
		*(index+i) = i;
	}
	omp_set_nested(1);
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
				//printf("Inner %d - %f\n",in+innersize*j,*(inner+in+innersize*j) );

			}
			in++;
		}
		else
		{
			*(outerindex+out) = *(index+i);
			for(int j=0;j<d;j++)
			{
				*(outer+out+outersize*j) = *(X+i+n*j);
				//printf("Outer %d - %f\n", out+outersize*j, *(outer+out+outersize*j));
			}
			out++;
		}
	}


	//Print Nodes of VP Tree
/*	printf("\t\t=======%d=======\n",n );
	printf("\t\tidx: %d\n",T->idx );
	printf("\t\tmedian: %f\n", T->md );
	printf("\t\tvp: %f %f \n",T->vp[0],T->vp[1]);*/
	if(omp_get_active_level()<5)
	{	
		#pragma omp parallel
			#pragma omp single
		{
			vptree *temp1,*temp2;
			temp1 = malloc(sizeof(vptree));
			temp2 = malloc(sizeof(vptree));

			#pragma omp task shared(temp1)
			temp1 = recbuildvp(inner,innerindex,innersize,d);
			#pragma omp task shared(temp2)
			temp2 = recbuildvp(outer,outerindex,outersize,d);
			#pragma omp taskwait
			T->inner = temp1;
			T->outer = temp2;
		}
		 
	}
	else
	{
		T->inner = recbuildvp(inner,innerindex,innersize,d);
		T->outer = recbuildvp(outer,outerindex,outersize,d);
	}
	
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
		#pragma omp parallel for schedule(static)
		for(int i=0;i<n-1;i++)
		{
			for(int j=0;j<d;j++)
			{
				*(distance+i) += (*(X+n*j+i) - *(X+n*j+n-1))*(*(X+n*j+i) - *(X+n*j+n-1));
			}
			*(distance+i) = sqrt(*(distance+i));
			//printf("%f\n", distance[i]);
		}
		return distance;
	}
	for(int i=0;i<n-1;i++)
	{
		for(int j=0;j<d;j++)
		{
			*(distance+i) += (*(X+n*j+i) - *(X+n*j+n-1))*(*(X+n*j+i) - *(X+n*j+n-1));
		}
		*(distance+i) = sqrt(*(distance+i));
		//printf("%f\n", distance[i]);
	}
	return distance;

}

double median(double *distance, int n)
{
	double md;
	if((n-1)%2!=0)
	{
		md = ksmallest(distance,n-1,(n-1)/2+1);
		//printf("\nMedian = %f\n", md);
	}
	else
	{
		md = (ksmallest(distance,n-1,(n-1)/2+1)+ksmallest(distance,n-1,(n-1)/2))/2;
		//printf("\nMedian = %f\n", md);
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