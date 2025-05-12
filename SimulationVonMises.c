/*
*	Jeff Peitsch
*
*	May 2025
*
*	Simulates epidemic using directional ILM with von Mises distribution
* 	on a population of 200 individuals for 20 time points
*
*/

#include <stdio.h> //for file io
#include <stdlib.h>
#include <time.h> //for random number
#include <math.h>

#define N 200 //size of population


double pdf(double phi, float kappa)//Von Mises distribution 
{
	float theta = 0.0;
	
	float t = kappa/3.75; //only accurate for kappa <= 3.75
	
	double bessel = 1.0 + 3.5156229*t*t + 3.0899424*pow(t,4.0) + 1.2067492*pow(t,6.0)
					 + 0.2659732*pow(t,8.0) + 0.0360768*pow(t,10.0) + 0.0045813*pow(t,12.0);
	
	return 1.0/2.0/M_PI/bessel*exp(kappa*cos(phi-theta));
}


void epigen(double x[N], double y[N], double alpha, double beta, double epsilon, double gamma, double kappa, int infections[N])
{
	double theta;
	
	for(int i = 0; i < N; i++)
		infections[i] = 0;
		
	infections[190] = 1;//set initial infection 
	
	
	for(int t = 1; t <= 20; t++)//loop through all time points
	{
		for(int i = 0; i < N; i++)//loop through all susceptible individuals
		{
			if(infections[i] == 0)
			{
			 	double dij = 0.0;
			 	for(int j = 0; j < N; j++)//loop through all infectious individuals
			 	{
			 		if(infections[j] != 0  && (infections[j] <= t  && infections[j]+gamma > t))
			 		{
			 			theta = atan2((y[i]-y[j]),(x[i]-x[j]));
			 			dij =  dij + pow(sqrt((x[i]-x[j])*(x[i]-x[j]) +(y[i]-y[j])*(y[i]-y[j])),(-beta/pdf(theta, kappa)));//infection kernel
			 		}
			 	}
			 	double prob = 1 - exp(-(alpha*dij)-epsilon);//ILM equation (P(it))
			 	double u = (double)rand()/RAND_MAX;
			 	if(prob >= u)//decision to infect or not
			 		infections[i] = t+1;
			}
		}
	}
}




int main(void)
{
	double alpha = 0.1;
	double beta = 0.2;
	int gamma = 2;
	double epsilon = 0.0;
	double kappa = 0.5;
	
	
	int inftime[N];
	
	double x[N];
	double y[N];
	FILE *pop = fopen("population.txt", "r");//read in population x and y locations
	for(int i = 0; i < N; i++)
	{
		fscanf(pop, "%lf %lf\n", &x[i], &y[i]);
	}
	fclose(pop);
	
	
	srand(time(NULL));//random seed for random number
		
	char filename[100];
	for(int j=0; j<20; j++)
	{
		sprintf(filename, "train%d.txt", j);
	
		FILE *train = fopen(filename, "w");
		epigen(x, y, alpha, beta, epsilon, gamma, kappa, inftime);
		
		
		for(int t = 0; t < N; t++)
		{
			fprintf(train, "%f %f %d\n", x[t], y[t], inftime[t]);
		}
		fprintf(train, "\n");
	
		fclose(train);
	}
	
}








