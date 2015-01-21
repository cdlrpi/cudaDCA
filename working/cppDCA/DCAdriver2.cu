//
//	DCA Driver
//
//This file invokes all of the necessary function calls to prepare
//and simulate a compound pendulum system through the use of the
//recursive DCA algorithm.  Output is created in a format that is
//readable in python for answer checking and graphing purposes.

//Included Files
#include <malloc.h>
#include <iostream>
#include <math.h>
#include "funct_bin/classes.h"
#include "funct_bin/npy.h"
#include <math.h>
#include <fstream>
#include <limits>
#include <cuda.h>
//Function Prototypes
//	Function found in RK45.cu
void RK_45(double state[], double step, int n, InitBody *bs, Joint *js,double Y[]);

//	Functions found in Functs.cu
void pend_init(InitBody *bs,int n,double mass, double length);
void horizontal_drop(double x[],int n);
void set_up(double A[], double B[], double C[], int n , double h);

//Main function
int main()
{
	int	n=0;
	std::ofstream timedata;
	std::ofstream numbods;
	numbods.open("numbods.mtx");
	timedata.open("graph_cpu.mtx");
while(n<40000)
	{
		if(n<500)
		{
			n+=10;
		}
		else if( n<2000)
		{
			n+=100;
		}
		else if(n< 10000)
		{
			n+= 1000;
		}
		else if(n<40000)
		{
			n+=10000;
		}
	typedef std::numeric_limits< double > dbl;

	
	std::cout.precision(dbl::digits10);
	//myfile2<<tstep<<"  ";
	//Vals[0][0]=tstep;	//Save the length of a timestep for plotting purposes	
	//Vals[1][0]=tfinal;	//Save the final time for plotting purposes
	cudaEvent_t beginEvent;
	cudaEvent_t endEvent;
	cudaEventCreate( &beginEvent );
	cudaEventCreate( &endEvent );
	//Variable Declarations
	InitBody* bodies; //List of bodies used for initialization only
	Joint* joints;  //List of joints between bodies NOTE: This joint list is not used in this version
	double *inits;	//Initial conditions
	double *Y;	//Solution to each timestep
	
	//System Setup
	//int n=50; //Number of bodies
	bodies = new InitBody[n]; //List of initialization bodies is length n
	joints = new Joint[n]; //List of joints is length n
	inits = new double[2*n];	//Initial conditions are length 2*n
	Y = new double[2*n];	//Timestep solution is length 2*n
	pend_init(bodies,n,1.0,1.0); //Initialize mass, length, and inertia of all bodies 

	//Time Setup
	double tstep= 0.001; //Length of a timestep [s]
	double tfinal = 0.005; //Final time [s]
	int tlen = (int) floor(tfinal/tstep)+1;	//Number of timesteps

	//Matrix Output Setup
	int shape1[2] = { tlen , 2*n }, fortran_order = 0;	//Shape of solution matrix
	int shape2[2] = { 2 , n+1 };	//Shape of matrix holding information to calculate the energy
	

	//System Initialization
	horizontal_drop(inits,n);	//Set the initial conditions
	
	//Save the initial conditions in the solution matrix
	cudaEventRecord( beginEvent, 0 );
	//Numerical Integration
	for(int t=1; t<tlen; t++)	//Loop through every timestep
	{
		RK_45(inits,tstep,n,bodies,joints,Y);	//Find the solution at that timestep
		for(int i = 0; i<2*n;i++)	//Loop through the solution
		{
			inits[i]=Y[i];	//Use the solution as the initial conditions for the next timestep
		}
		
	}
	cudaEventRecord( endEvent, 0 );
	cudaEventSynchronize( endEvent );
	float timeValue;
	cudaEventElapsedTime( &timeValue, beginEvent, endEvent );
	timedata<< timeValue << "  ";
	numbods<<n<<"  ";
	std::cout << n << std::endl;
	
	
	//Free memory
	delete[] inits;
	delete[] Y;
	delete[] bodies;
	delete[] joints;
	}
	numbods.close();
	timedata.close();
	return EXIT_SUCCESS;	//Program completed successfully
}
