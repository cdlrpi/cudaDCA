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
#include "d_code/deviceDisassemble.h"
#include "d_code/deviceAssemble.h"
#include "d_code/deviceInitialize.h"
#include "d_code/deviceFuncts.h"
#include "funct_bin/npy.h"
#include <fstream>
#include <limits>
//Function Prototypes
//	Function found in RK45.cu
void RK_45(double state[], double step, int n, double m[], double l[], double II[],double Y[],int cut_off);

//	Functions found in Functs.cu
void pend_init(double m[], double l[], double II[],int n,double mass, double length);
void horizontal_drop(double x[],int n);
void set_up(double A[], double B[], double C[], int n , double h);

//Main function
int main()
{
	//Variable Declarations
	double *inits;	//Initial conditions
	double *Y;	//Solution to each timestep
	double *m;
	double *l;
	double *II;
	
	std::ofstream myfile;
	std::ofstream myfile2;
	myfile2.open("Vals.mtx");
  	myfile.open ("output.mtx");

	//System Setup
	int n=80000; //Number of bodies
	inits = (double*)malloc(sizeof(double)*2*n);	//Initial conditions are length 2*n
	Y = (double*)malloc(sizeof(double)*2*n);	//Timestep solution is length 2*n
	m = (double*)malloc(sizeof(double)*n);
	l = (double*)malloc(sizeof(double)*n);
	II = (double*)malloc(sizeof(double)*n*3*3);
	pend_init(m,l,II,n,1.0,1.0); //Initialize mass, length, and inertia of all bodies 

	/////////////////////////////////////////////////////
	int cut_off = 50;


	//Time Setup
	double tstep= 0.001; //Length of a timestep [s]
	double tfinal = 0.002; //Final time [s]
	int tlen = (int) floor(tfinal/tstep)+1;	//Number of timesteps
	
	typedef std::numeric_limits< double > dbl;

	std::cout.precision(15);
	myfile2<<tstep<<"  ";
	//Matrix Output Setup
	

	for(int i =1; i<n+1; i++)	//Loop through all of the bodies
	{
		myfile2<<m[i-1]<<"  ";
	}
	myfile2<<"\n"<<tfinal<<"  ";
	for(int i =1; i<n+1; i++)	//Loop through all of the bodies
	{
		myfile2<<l[i-1]<<"  ";
	}
	//System Initialization
	horizontal_drop(inits,n);	//Set the initial conditions
	
	//Save the initial conditions in the solution matrix
	for(int r=0; r<2*n; r++)
	{
		myfile << inits[r]<< "  ";
	}
	myfile << "\n";

	//Numerical Integration
	for(int t=1; t<tlen; t++)	//Loop through every timestep
	{
		RK_45(inits,tstep,n,m,l,II,Y,cut_off);	//Find the solution at that timestep
		//std::cout<<"hey"<<std::endl;
		for(int i = 0; i<2*n;i++)	//Loop through the solution
		{
			
			inits[i]=Y[i];	//Use the solution as the initial conditions for the next timestep
			myfile << inits[i]<<"  ";
			//std::cout<<"hey"<<std::endl;
		}
		myfile << "\n";
		
	}
	
	//Free memory
	free(inits);
	free(II);
	free(m);
	free(l);
	free(Y);
	myfile.close();
	myfile2.close();
	return EXIT_SUCCESS;	//Program completed successfully
}
