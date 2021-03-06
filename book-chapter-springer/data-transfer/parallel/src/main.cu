//
//	CUDA DCA Driver
//
//This file invokes all of the necessary function calls to prepare
//and simulate a compound pendulum system through the use of the
//recursive DCA algorithm.	The majority of this algorithm is run
//on the gpu.  Output is created in a format that is
//readable in python for answer checking and graphing purposes.

//Included Files
#include <malloc.h>
#include <iostream>
#include <math.h>
#include "dca-functs/classes.h"
#include "d-code/deviceDisassemble.h"
#include "d-code/deviceAssemble.h"
#include "d-code/deviceInitialize.h"
#include "d-code/deviceFuncts.h"
#include <math.h>
#include <fstream>
#include <limits>


//Function Prototypes
//	Function found in RK45.cu
void rk45(double state[], double step, int n, double m[], double l[], double II[],double Y[],int cut_off, double Zs[], float times[], int reps);

//	Functions found in Functs.cu
void pend_init(double m[], double l[], double II[],int n,double mass, double length);
void horizontal_drop(double x[],int n);
void set_up(double A[], double B[], double C[], int n , double h);
void Initialize(double m[], double l[],double II[],double Zetas[],int n);

//Main function
int main()
{
	//reps is the number of times DCA will repeat itself at the current
	//number of bodies.  If reps is set to 5, DCA will be timed 5 times
	//and the average of the 5 runs is the output.  The output time will
	//always be the time in ms for a single DCA run.
	int reps=1;
	

	int	n=0;
	int cut_off;
	float *times=(float*)malloc(sizeof(float)*reps);
	std::ofstream timedata;
	std::ofstream numbods;

//FILE NAMES
//numbods is a list of the number of bodies used for each run
//timedata is a matrix that holds the time it took for each run
	numbods.open("n.mtx");
	timedata.open("t.mtx");

////////////////////////////////////////////////////////////////
	int numa;
	std::cout<<"\nStarting Speed Test..."<<std::endl;
//This loop determines the number of assemblies (numa) to do on the gpu
//Right now it is set to do 5 runs total, the first run does no assemblies
//on the gpu, the second does 1 assembly, the third does 3, the fourth does 6, the fith does all levels on the GPU.
//You can change these numbers however you want and the code will adapt and only 
//do as many as is needed (if you ask for 12 assemblies on 2 bodies it will still
//only assemble once)
	for(int xx = 0; xx<5; xx+=1)
	{
		if(xx ==0)//This should have been a switch statement
		{
			numa = 0;
		}
		if(xx == 1)
		{
			numa = 1;
		}
		if (xx ==2)
		{
			numa = 3;
		}
		if(xx ==3)
		{
			numa = 6;
		}
		if(xx ==4)
		{
			numa = 12;
			// for 2048 bodies this should be all levels
		}
	n = 2;

//This loop cycles from 0 to the desired maximum number of bodies (in this case 2048)
//The if statements inside determine how the number of bodies, n , is incremented.
//This is because you may need points that are closer together for the spots in the graph where 
//curves happen.  You can change the increment however you want to get the point density you need
//You can also put a "reps =" somthing line in each if statement if you want to take a bigger
//average of points at first, but don't feel like waiting for larger numbers of bodies.
	while(n<2048)
	{
		if(n<512)
		{
			n+=1;
			//reps = 20;
		}
		else if( n<2048)
		{
			n+=64;
			//reps = 10;
		}
		else if(n<8192)
		{
			n+= 512;
			//reps = 5;
		}
		else
		{
			n+=2048;
		}

		float *times=(float*)malloc(sizeof(float)*reps);
		int x = n;
		times[0] = 0;
		cut_off=x;
		for(int c =0; c<numa; c++)
		{
			if(x==1)
			{
			cut_off=0;
			
		}
		else if( x%2==0)
		{
			x=x/2;
		}
		else
		{
			x++;
			x=x/2;
		}
		}
	

	if(cut_off !=0)
	{
		cut_off =x;
	}
	
	//Variable Declarations
	//Variable Declarations
	double *inits;	//Initial conditions
	double *Y;	//Solution to each timestep
	double *m;
	double *l;
	double *II;  
	double *Zs;

	//System Setup
	inits = (double*)malloc(sizeof(double)*2*n);	//Initial conditions are length 2*n
	Y = (double*)malloc(sizeof(double)*2*n);	//Timestep solution is length 2*n
	m = (double*)malloc(sizeof(double)*n);
	l = (double*)malloc(sizeof(double)*n);
	II = (double*)malloc(sizeof(double)*n*3*3);
	Zs = (double*)malloc(sizeof(double)*n*26*6);
	pend_init(m,l,II,n,1.0,1.0); //Initialize mass, length, and inertia of all bodies 

	Initialize(m, l, II, Zs, n); //Find initial (unrotated) zeta values

	/////////////////////////////////////////////////////

	//Time Setup
	double tstep= 0.001; //Length of a timestep [s]
	double tfinal = 0.001; //Final time [s] set to only perform one DCA run for timing purposes
	int tlen = (int) floor(tfinal/tstep)+1;	//Number of timesteps	

	//Set precision of output
	typedef std::numeric_limits< double > dbl;
	std::cout.precision(dbl::digits10);

	//System Initialization
	horizontal_drop(inits,n);	//Set the initial conditions

	rk45(inits,tstep,n,m,l,II,Y,cut_off,Zs,times, reps);//Call the integrator, in this version the integrator is not functioning for timing purposes

	//find the average of the recorded times
	float timeValue=0;
	for(int i =0; i<reps; i++)
	{
		timeValue+=times[i];
	}
	timeValue/=reps;
	timedata<< timeValue << "  ";
	numbods<<n<<"  ";

	//Check for an error, if one occurs there is a problem with your CUDA installation
 	if ( cudaSuccess != cudaGetLastError() ){
	std::cerr<<"Error!  Problem with CUDA installation!"<<std::endl;
	exit(1);}
	
	//Free memory
	free(inits);
	free(II);
	free(m);
	free(l);
	free(Y);
	}
	std::cout<<"Serial Test Complete."<<std::endl;
	timedata<<"\n";
	numbods<<"\n";

	}
	
	numbods.close();
	timedata.close();
	return EXIT_SUCCESS;	//Program completed successfully
}	
