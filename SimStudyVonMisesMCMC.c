/*
*	Jeff Peitsch
*
*	May 2025
*
*	Calculation of Bayesian posterior distribution 
* 	for the directional ILM with von Mises distribution
*	using random walk MH-MCMC
*
*/

#include <stdio.h> //for file io
#include <stdlib.h>
#include <time.h> //for random number and run timer
#include <math.h>
#include <float.h>



#define LENGTH 20 //length of epidemic
#define N 200 // number of individuals

#define ITER 20000


//population data variables:
int infections[N];
float x[N];
float y[N];


//prior distribution functions:
double alphaprior(float z);
double thetaprior(float z);
double kappaprior(float z);
double betaprior(float z);

//log-likelihood function
double loglike(float alpha, float beta, float theta, float kappa);


int main(void)
{
	float u,u1,u2,u3,u4,psi;
	float aout[ITER];
	float bout[ITER];
	float tout[ITER];
	float kout[ITER];
	double like[ITER];
	
	
	srand(time(NULL));//random seed for random number
	
	
	//!!!!!!!!Proposal parameter !!!!!!!!!!!!!!!!!!!!!!!
	float aa=0.02;
	float bb=0.02;
	float tt=0.8;
	float kk=0.1;
	
	
	
	//repeat entire program for 20 different epidemics
	char filename[100];
for(int file = 0; file <20; file++)
{
	sprintf(filename, "train%d.txt", file);
	FILE *pop = fopen(filename, "r");
	for(int i = 0; i < N; i++)//read epidemic dataset
	{
		fscanf(pop, "%f %f %d", &x[i], &y[i], &infections[i]);
	}
	fclose(pop);
	
	
	

	//starting values 
	aout[0]=0.1;
	bout[0] = 0.2;
	tout[0]=0.0;
	kout[0] = 0.3;

	
	
	
	int rejectabk = 0;
	int rejectt = 0;
	
	
	//!!!!!!!!!!!!!!!!!!MCMC Chain!!!!!!!!!!!!!!!!


	//! Start timer
	time_t begin = time(NULL);
	
	double liketemp;


	like[0] = loglike(aout[0], bout[0], tout[0], kout[0]);
	

	for(int i = 0; i < ITER-1; i++)
	{
		
		//rand() generates a number between 0 and RAND_MAX
		u1 = 2*aa*((float)rand()/RAND_MAX) - aa;
		aout[i+1] = aout[i]+u1;
		
		u2 = 2*bb*((float)rand()/RAND_MAX) - bb;
		bout[i+1] = bout[i]+u2;
		
		u3 = 2*tt*((float)rand()/RAND_MAX) - tt;
		tout[i+1] = tout[i]+u3;
		
		u4 = 2*kk*((float)rand()/RAND_MAX) - kk;
		kout[i+1] = kout[i]+u4;
		
		
   		
		//Update alpha, beta, kappa in a block
		 if(aout[i+1] <= 0.0 || bout[i+1] <= 0.0 || kout[i+1] < 0.0) 
	    {
	    	aout[i+1]= aout[i];
		  	like[i+1]=like[i];
		  	bout[i+1] = bout[i];
		  	kout[i+1] = kout[i];
		  	rejectabk += 1;
		}
	    else
		{		
			
			like[i+1] = loglike(aout[i+1], bout[i+1], tout[i], kout[i+1]);
			
		   	psi = (like[i+1]-like[i])+(alphaprior(aout[i+1])-alphaprior(aout[i]))
		   							+(betaprior(bout[i+1])-betaprior(bout[i]))+ 
		   							(kappaprior(kout[i+1])-kappaprior(kout[i]));
			psi = 0.0<psi? 0.0 : psi; // minimum function
			
				
			u=log((float)rand()/RAND_MAX);
				
		  	//!!!!!!!!!!!!!! Test for acceptance probability!!!!!!!!!!!!!!!!!!!!! 
	  
			if(psi < u)
			{
				aout[i+1]= aout[i];
				like[i+1]=like[i];
				bout[i+1]= bout[i];
		  		kout[i+1]= kout[i];
		  	rejectabk += 1;
			}
	   	}

		liketemp = like[i+1];
		
		
		//Update theta in its own
	    if(0) 
	    {
		  	tout[i+1]= tout[i];
		  	rejectt += 1;
		}
	    else
		{		
			if(tout[i+1] < -M_PI)
			{
				tout[i+1] += 2*M_PI;
			}
			if(tout[i+1] > M_PI)
			{
				tout[i+1] -= 2*M_PI;
			}
			
			like[i+1] = loglike(aout[i+1], bout[i+1], tout[i+1], kout[i+1]);
			
		   	psi = (like[i+1]-liketemp)+ (thetaprior(tout[i+1])-thetaprior(tout[i]));
			psi = 0.0<psi? 0.0 : psi; // minimum function
			
				
			u=log((float)rand()/RAND_MAX);
				
		  	//!!!!!!!!!!!!!! Test for acceptance probability!!!!!!!!!!!!!!!!!!!!! 
	  
			if(psi < u)
			{
				tout[i+1]= tout[i];
				like[i+1]=liketemp;
		  	rejectt += 1;
			}
	   	}



	}

	//Output MCMC chain to file
	sprintf(filename, "output%d.csv", file);
	
	FILE *out = fopen(filename, "w");
	fprintf(out, "alpha,beta,theta,kappa,log-like\n");
	for(int i = 0; i < ITER; i++)
	{
		fprintf(out, "%f,%f,%f,%f,%lf\n", aout[i], bout[i], tout[i], kout[i], like[i]);
	}
	fclose(out);
 
 	//print rejection rates
	printf("abk: %f\n", (float)rejectabk/ITER);
	printf("theta: %f\n", (float)rejectt/ITER);
 
// 	! End timer
	time_t end = time(NULL);

	//Calculate run-time
	float tdiff_s = (float)(end - begin)/60;
	
	printf("That took %f minutes\n", tdiff_s);

}
}



//!!!!!!!!!!!!!Prior!!!!!!!!!!!!!!!!!!!!!!!!

///////////all priors return the log value/////////////


double betaprior(float z)//exponential(mean = 10000)
{
  	if(z < 0.0)
     	return -FLT_MAX;
	
	float lambda = 10000;
	
	return -z/lambda;
}

double alphaprior(float z)//exponential(mean = 10000)
{
  	if(z < 0.0)
     	return -FLT_MAX;
     	
	float lambda = 10000;
	
	return -z/lambda;
}


double thetaprior(float z)//uniform(-pi, pi)
{
  	if(z < -M_PI || z > M_PI)
     	return -FLT_MAX;
  	 
  	 return 0;//only need up to proportionality
}


//for Von Mises distribution
double kappaprior(float z)//exponential(mean = 10000)
{
  	if(z < 0.0)
     	return -FLT_MAX;  	 
	
	float lambda = 10000;
	
	return -z/lambda;

}


double pdf(double phi, float theta, float kappa)//Von Mises distribution 
{
	float t = kappa/3.75; //only accurate for kappa <= 3.75
	
	double bessel = 1.0 + 3.5156229*t*t + 3.0899424*pow(t,4.0) + 1.2067492*pow(t,6.0)
					 + 0.2659732*pow(t,8.0) + 0.0360768*pow(t,10.0) + 0.0045813*pow(t,12.0);
	
	return 1.0/2.0/M_PI/bessel*exp(kappa*cos(phi-theta));
}



double loglike(float alpha, float beta, float theta, float kappa)
{	
	//parameters being held constant:
	float epsilon = 0.0;
	int gamma = 2;
	
	if(alpha < 0 || beta < 0 || kappa < 0)
	{
		printf("flag");
	    return -FLT_MAX;//0.0;
	}
	
	if(theta < -M_PI)
	{
		theta += 2.0*M_PI;
	}
	else if(theta > M_PI)
	{
		theta -= 2.0*M_PI;
	}
	
    double like = 0.0;

	
	for(int t = 1; t < LENGTH; t++)//loop through all time points
	{
		for(int i = 0; i < N; i++)//loop through all susceptible individuals
		{
			if(infections[i] == 0 || infections[i] > t)
			{
			 	long double dij = 0.0;
			 	for(int j = 0; j < N; j++)//loop through all infectious individuals
			 	{
			 		if(infections[j] != 0  && (infections[j] <= t  && infections[j]+gamma > t))
			 		{
			 			double phi = atan2((y[i]-y[j]),(x[i]-x[j]));
			 			dij += pow(sqrt((x[i]-x[j])*(x[i]-x[j]) +(y[i]-y[j])*(y[i]-y[j])),(-beta/pdf(phi, theta, kappa)));//infection kernel
			 		}
			 	}
			 	long double prob = (long double)1.0 - exp(-((long double)alpha*dij)-epsilon);//ILM equation
			 	
				if(infections[i] == 0 || infections[i] > t+1)//susceptible individuals
				{
					like += log(1.0-prob);
				}
				else if(infections[i]==t+1)//newly infected individuals
				{
					like += log(prob);
				}
			}
		}
	}
    
    return like;
}
