/*
*	Jeff Peitsch
*
*	May 2025
*
*	Calculation of Bayesian posterior distribution 
* 	for the directional FMD-ILM with wrapped Cauchy distribution
*	using random walk MH-MCMC
*
*/

#include <stdio.h> //for file io
#include <stdlib.h>
#include <time.h> //for random number and run timer
#include <math.h>
#include <float.h>



#define LENGTH 60 //length of epidemic
#define N 1419 // number of individuals

#define ITER 80000


//population data variables:
int infections[N];
int x[N];
int y[N];
int cattle[N];
int sheep[N];


//prior distribution functions:
double alphasprior(float z);
double alphacprior(float z);
double thetaprior(float z);
double kappaprior(float z);
double betaprior(float z);
double epsilonprior(float z);

//log-likelihood function
double loglike(float alpha_s,float alpha_c, float beta, float theta, float kappa, float epsilon);


int main(void)
{
	float u,u1,u2,u3,u4,psi;
	float asout[ITER];
	float acout[ITER];
	float bout[ITER];
	float tout[ITER];
	float kout[ITER];
	float eout[ITER];
	double like[ITER];
	
	
	srand(time(NULL));//random seed for random number
	
	
	//!!!!!!!!Proposal parameter !!!!!!!!!!!!!!!!!!!!!!!
	float as=0.005;
	float ac=0.06;
	float bb=0.01;
	float tt=0.9;
	float kk=0.06;
	float ee = 0.0005;
	
	
	
	int county;
	int p;
	int h;
	int region;
	int area;
	int pigs;
	int goats;
	int other;
	
	
	FILE *pop = fopen("FMDpop.txt", "r");
	for(int i = 0; i < N; i++)
	{
		fscanf(pop, "%d%d%d%d%d%d%d%d%d%d%d%d\n", &county, &p, &h, &x[i], &y[i], &area, 
				&cattle[i], &pigs, &sheep[i], &goats, &other, &infections[i]);
	}
	fclose(pop);
	
	

	//initial values
	asout[0]=0.001;
	acout[0]=0.1;
	bout[0] = 0.2;
	tout[0]=3.0;
	kout[0] = 0.05;
	eout[0] = 0.0001;

	
	
	
	int rejectas = 0;
	int rejectac = 0;
	int rejectb = 0;
	int rejectk = 0;
	int rejectt = 0;
	int rejecte = 0;
	
	
	//!!!!!!!!!!!!!!!!!!MCMC Chain!!!!!!!!!!!!!!!!


	//! Start timer
	time_t begin = time(NULL);
	
	double liketemp;


	like[0] = loglike(asout[0], acout[0], bout[0], tout[0], kout[0], eout[0]);
	

	for(int i = 0; i < ITER-1; i++)
	{
		
		//rand() generates a number between 0 and RAND_MAX
		u1 = 2*as*((float)rand()/RAND_MAX) - as;
		asout[i+1] = asout[i]+u1;
		
		u1 = 2*ac*((float)rand()/RAND_MAX) - ac;
		acout[i+1] = acout[i]+u1;
		
		u2 = 2*bb*((float)rand()/RAND_MAX) - bb;
		bout[i+1] = bout[i]+u2;
		
		u3 = 2*tt*((float)rand()/RAND_MAX) - tt;
		tout[i+1] = tout[i]+u3;
		
		u4 = 2*kk*((float)rand()/RAND_MAX) - kk;
		kout[i+1] = kout[i]+u4;
		
		u4 = 2*ee*((float)rand()/RAND_MAX) - ee;
		eout[i+1] = eout[i]+u4;
		
   		
		
		 if(asout[i+1] <= 0.0) 
	    {
	    	asout[i+1]= asout[i];
		  	like[i+1]=like[i];
		  	rejectas += 1;
		}
	    else
		{		
			like[i+1] = loglike(asout[i+1], acout[i], bout[i], tout[i], kout[i], eout[i]);
			
		   	psi = (like[i+1]-like[i])+(alphasprior(asout[i+1])-alphasprior(asout[i]));
			psi = 0.0<psi? 0.0 : psi; // minimum function
			
			u=log((float)rand()/RAND_MAX);
				
		  	//!!!!!!!!!!!!!! Test for acceptance probability!!!!!!!!!!!!!!!!!!!!! 
			if(psi < u)
			{
				asout[i+1]= asout[i];
				like[i+1]=like[i];
		  	rejectas += 1;
			}
	   	}
	   	
		liketemp = like[i+1];
	   	
	   	if(acout[i+1] <= 0.0 || bout[i+1] <= 0.0) 
	    {
	    	acout[i+1]= acout[i];
		  	like[i+1]=liketemp;
		  	bout[i+1] = bout[i];
		  	rejectac += 1;
		}
	    else
		{		
			like[i+1] = loglike(asout[i+1], acout[i+1], bout[i+1], tout[i], kout[i], eout[i]);
			
		   	psi = (like[i+1]-liketemp)+(alphacprior(acout[i+1])-alphacprior(acout[i]))+
		   								(betaprior(bout[i+1])-betaprior(bout[i]));
			psi = 0.0<psi? 0.0 : psi; // minimum function
			
			u=log((float)rand()/RAND_MAX);
				
		  	//!!!!!!!!!!!!!! Test for acceptance probability!!!!!!!!!!!!!!!!!!!!! 
			if(psi < u)
			{
				acout[i+1]= acout[i];
				bout[i+1]= bout[i];
				like[i+1]=liketemp;
		  	rejectac += 1;
			}
	   	}
	   	
		liketemp = like[i+1];
	   
	   
	   	if(kout[i+1] < 0.0) 
	    {
		  	like[i+1]=liketemp;
		  	kout[i+1] = kout[i];
		  	rejectk += 1;
		}
	    else
		{		
			
			like[i+1] = loglike(asout[i+1], acout[i+1], bout[i+1], tout[i], kout[i+1], eout[i]);
			
		   	psi = (like[i+1]-liketemp)+(kappaprior(kout[i+1])-kappaprior(kout[i]));
			psi = 0.0<psi? 0.0 : psi; // minimum function
			
				
			u=log((float)rand()/RAND_MAX);
				
		  	//!!!!!!!!!!!!!! Test for acceptance probability!!!!!!!!!!!!!!!!!!!!! 
	  
			if(psi < u)
			{
				like[i+1]=liketemp;
		  		kout[i+1]= kout[i];
		  	rejectk += 1;
			}
	   	}

		liketemp = like[i+1];
		
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
			
			like[i+1] = loglike(asout[i+1], acout[i+1], bout[i+1],tout[i+1], kout[i+1], eout[i]);
			
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
	   	
		liketemp = like[i+1];
	   	
	   	if(eout[i+1] < 0.0) 
	    {
		  	like[i+1]=liketemp;
		  	eout[i+1] = eout[i];
		  	rejecte += 1;
		}
	    else
		{		
			
			like[i+1] = loglike(asout[i+1], acout[i+1], bout[i+1], tout[i+1], kout[i+1], eout[i+1]);
			
		   	psi = (like[i+1]-liketemp)+(epsilonprior(eout[i+1])-epsilonprior(eout[i]));
			psi = 0.0<psi? 0.0 : psi; // minimum function
			
				
			u=log((float)rand()/RAND_MAX);
				
		  	//!!!!!!!!!!!!!! Test for acceptance probability!!!!!!!!!!!!!!!!!!!!! 
	  
			if(psi < u)
			{
				like[i+1]=liketemp;
		  		eout[i+1]= eout[i];
		  	rejecte += 1;
			}
	   	}
	}

	
	FILE *out = fopen("output.csv", "w");
	fprintf(out, "alpha_s,alpha_c,beta,theta,kappa,epsilon,log-like\n");
	for(int i = 0; i < ITER; i++)
	{
		fprintf(out, "%f,%f,%f,%f,%f,%f,%lf\n", asout[i], acout[i], bout[i], tout[i], kout[i], eout[i], like[i]);
	}
	fclose(out);
 
	printf("as: %f\n", (float)rejectas/ITER); 
	printf("ac: %f\n", (float)rejectac/ITER); 
	printf("b: %f\n", (float)rejectb/ITER); 
	printf("k: %f\n", (float)rejectk/ITER);
	printf("e: %f\n", (float)rejecte/ITER);
	printf("theta: %f\n", (float)rejectt/ITER);
 
// 	! End timer
	time_t end = time(NULL);

	//Calculate run-time
	double tdiff_s = (double)(end - begin)/60.0/60.0;
	
	printf("That took %lf hours\n", tdiff_s);


}



//!!!!!!!!!!!!!Prior!!!!!!!!!!!!!!!!!!!!!!!!

///////////all priors return the log value/////////////


double betaprior(float z)//exponential(1)
{
  	if(z < 0.0)
     	return -FLT_MAX;
	
	float lambda = 1;
	
	return -z/lambda;
}

double alphasprior(float z)//exponential(100)
{
  	if(z < 0.0)
     	return -FLT_MAX;
     	
	float lambda = 0.01;
	
	return -z/lambda;
}

double alphacprior(float z)//exponential(1)
{
  	if(z < 0.0)
     	return -FLT_MAX;
     	
	float lambda = 1;
	
	return -z/lambda;
}


double thetaprior(float z)//uniform(-pi, pi)
{
  	if(z < -M_PI || z > M_PI)
     	return -FLT_MAX;
  	 
  	 return 0;//only need up to proportionality
}


//for wrapped Cauchy distribution
double kappaprior(float z)//uniform(0,1)
{
  	if(z < 0.0 || z > 1.0)
     	return -FLT_MAX;//0.0;
     	
    return 0;
}


double epsilonprior(float z)//exponential(10)
{
  	if(z < 0.0)
     	return -FLT_MAX;
  	 
	
	float lambda = 0.1;
	
	return -z/lambda;

}




double pdf(double phi, float theta, float kappa)//Wrapped Cauchy distribution 
{
	
	return 1.0/2.0/M_PI*(1.0-kappa*kappa)/(1.0+kappa*kappa-2.0*kappa*cos(phi-theta));
}



double loglike(float alpha_s, float alpha_c, float beta, float theta, float kappa, float epsilon)
{	
	//parameters being held constant:
	int gamma = 4;
	
	if(alpha_c < 0 || beta < 0 || kappa < 0)
	{
	    return -FLT_MAX;
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
			 		if(infections[j] > 0  && (infections[j] <= t  && infections[j]+gamma > t))
			 		{
			 			double phi = atan2((y[i]-y[j]),(x[i]-x[j]));
			 			dij += pow(sqrt((double)(x[i]-x[j])*(x[i]-x[j]) +(double)(y[i]-y[j])*(y[i]-y[j])),(-beta/pdf(phi, theta, kappa)));//infection kernel
			 		}
			 	}
			 	long double prob = 1 - exp(-(long double)(alpha_s*sheep[i]+alpha_c*cattle[i])*dij-epsilon);//FMD-ILM equation
			 	
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
