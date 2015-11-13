//
//	Hybrid DCA Driver
//
//This file invokes all of the necessary function calls to prepare
//and simulate a compound pendulum system through the use of the
//recursive DCA algorithm.	Much of this algorithm is run
//on the gpu.

//Defined values to interperate flag
#define python_file 0
#define mtx_file 1
#define both 2


//Included Files
//	standard files
#include <malloc.h>
#include <iostream>
#include <string.h>
#include <math.h>
#include <cuda.h>
#include <fstream>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
//	created files
#include "Device_Code/deviceDisassemble.h"
#include "Device_Code/deviceAssemble.h"
#include "Device_Code/deviceInitialize.h"
#include "Device_Code/deviceFuncts.h"
#include "DCA_Functions/npy.h"

//Function Prototypes
//	Functions found in RK45.cu
void RK_45(double state[], double step, int n, double m[], double l[], double II[],double Y[],int cut_off, double Zs[]);

//	Functions found in Functs.cu
void pend_init(double m[], double l[], double II[],int n,double mass, double length);
void horizontal_drop(double x[],int n);
void set_up(double A[], double B[], double C[], int n , double h);

//	Functions found in Igaoibn;aeirn
void Initialize(double m[], double l[],double II[],double Zetas[],int n);

// Functions found in findCutoff.cu
int findCutoff(int n, int accuracy);

using namespace std;	//namespace for output and input

//Main function
int main(int argc, char* argv[])
{
	//Variable Declarations
	double *inits;	//Initial conditions
	double *Y;	//Solution to each timestep
	double *m;	//Mass of each body
	double *l;	//Length of each body
	double *II;	//Inertia matrix of each body
	double *Zs;	//Zeta matrices for all bodies
	double tfinal;	//Final time
	int n;	//Number of bodies
	int flag;	//Flag to keep track of output
	int numa;	//Number of assemblies on gpu
	int percent = 10; //Number for printing progress
	int tlen;	//Number of timesteps
	double percheck = 0.1;	//Number to check progress
	double tstep=.001;	//Length of a timestep
	char c;
	float timeValue;

	//reps is the number of times DCA will repeat itself at the current
	//number of bodies.  If reps is set to 5, DCA will be timed 5 times
	//and the average of the 5 runs is the output.  The output time will
	//always be the time in ms for a single DCA run.
	int reps=1;
	

	n=0;
	int cut_off;
	float *times=(float*)malloc(sizeof(float)*reps);
	std::ofstream timedata;
	std::ofstream numbods;

//FILE NAMES
//numbods is a list of the number of bodies used for each run
//timedata is a matrix that holds the time it took for each run
	numbods.open("serialBodies.mtx");
	timedata.open("serialTimes.mtx");

////////////////////////////////////////////////////////////////

//This loop determines the number of assemblies (numa) to do on the gpu
//Right now it is set to do 4 runs total, the first run does no assemblies
//on the gpu, the second does 1 assembly, the third does 3, and the fourth does 6.
//You can change these numbers however you want and the code will adapt and only 
//do as many as is needed (if you ask for 12 assemblies on 2 bodies it will still
//only assemble once)
	for(int xx = 0; xx<1; xx+=1)
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
		}
		
	n=0;
	std::cout<<"\n\n\n\n\n"<<numa<<"\n\n\n\n";

//This loop cycles from 0 to the desired maximum number of bodies (in this case 4000)
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
		else if(n< 8192)
		{
			n+= 512;
			//reps = 5;
		}
		else
		{
			n+=10000;
		}


	//Time Setup
	double tstep= 0.001; //Length of a timestep [s]
	double tfinal = 0.20; //Final time [s]
	int tlen = (int) floor(tfinal/tstep)+1;	//Number of timesteps


	//Allocate memory
	inits = (double*)malloc(sizeof(double)*2*n);	//Initial conditions are length 2*n
	Y = (double*)malloc(sizeof(double)*2*n);	//Timestep solution is length 2*n
	m = (double*)malloc(sizeof(double)*n);
	l = (double*)malloc(sizeof(double)*n);
	II = (double*)malloc(sizeof(double)*n*3*3);
	Zs = (double*)malloc(sizeof(double)*n*26*6);
	
	pend_init(m,l,II,n,1.0,1.0); //Initialize mass, length, and inertia of all bodies 
	Initialize(m, l, II, Zs, n); //Initialize the Zeta matrices
	horizontal_drop(inits,n);	//Set the initial conditions
	
	//Set maximum output precision
	std::cout.precision(16);
	std::cout<<"n = "<<n<<std::endl;
	cudaEvent_t beginEvent;
	cudaEvent_t endEvent;
	cudaEventCreate( &beginEvent );
	cudaEventCreate( &endEvent );

	cudaEventRecord( beginEvent, 0 );
	for(int t=1; t<tlen; t++)	//Loop through every timestep
	{
		RK_45(inits,tstep,n,m,l,II,Y,numa,Zs);	//Find the solution at that timestep
		for(int i = 0; i<2*n;i++)	//Loop through the solution
		{
			inits[i]=Y[i];	//Use the solution as the initial conditions for the next timestep
		}

	}
	cudaEventRecord( endEvent, 0 );
	cudaEventSynchronize( endEvent );
	cudaEventElapsedTime( &timeValue, beginEvent, endEvent );
	cout<<"100%\n\t\t";

	//Check for errors in integration and complete output
	if ( cudaSuccess != cudaGetLastError() )
	{
	 	cout<<"\nAn error occurred during integration, No output generated\n";
		return 0;
	}
	timedata<< timeValue << "  ";
	numbods<<n<<"  ";
	free(inits);	//Initial conditions
	free(Y);	//Solution to each timestep
	free(m);
	free(l);
	free(II);
	free(Zs);

		}
	timedata<<"\n";
	numbods<<"\n";
	}
	
	numbods.close();
	timedata.close();
	return EXIT_SUCCESS;	//Program completed successfully
}
