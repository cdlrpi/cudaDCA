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
#include <cuda.h>
#include "d_code/deviceDisassemble.h"
#include "d_code/deviceAssemble.h"
#include "d_code/deviceInitialize.h"
#include "d_code/deviceFuncts.h"
#include "funct_bin/npy.h"

//Function Prototypes
//	Function found in RK45.cu
void RK_45(double state[], double step, int n, double m[], double l[], double II[],double Y[],int cut_off, double Zs[]);

//	Functions found in Functs.cu
void pend_init(double m[], double l[], double II[],int n,double mass, double length);
void horizontal_drop(double x[],int n);
void set_up(double A[], double B[], double C[], int n , double h);
void Initialize(double m[], double l[],double II[],double Zetas[],int n);
//Main function
int main()
{
	//Variable Declarations
	double *inits;	//Initial conditions
	double *Y;	//Solution to each timestep
	double *m;
	double *l;
	double *II;
	double *Zs;
	
	//System Setup
	int n=167; //Number of bodies
	inits = (double*)malloc(sizeof(double)*2*n);	//Initial conditions are length 2*n
	Y = (double*)malloc(sizeof(double)*2*n);	//Timestep solution is length 2*n
	m = (double*)malloc(sizeof(double)*n);
	l = (double*)malloc(sizeof(double)*n);
	II = (double*)malloc(sizeof(double)*n*3*3);
	Zs = (double*)malloc(sizeof(double)*n*26*6);
	pend_init(m,l,II,n,1.0,1.0); //Initialize mass, length, and inertia of all bodies 

	Initialize(m, l, II, Zs, n);
	

	/////////////////////////////////////////////////////
	int cut_off = 50;


	//Time Setup
	double tstep= 0.001; //Length of a timestep [s]
	double tfinal =2; //Final time [s]
	int tlen = (int) floor(tfinal/tstep)+1;	//Number of timesteps

	//Matrix Output Setup
	int shape1[2] = { tlen , 2*n }, fortran_order = 0;	//Shape of solution matrix
	int shape2[2] = { 2 , n+1 };	//Shape of matrix holding information to calculate the energy
	double Vals[2][n+1];	//Matrix holding information to calculate the energy
	double X[tlen][2*n];	//Solution Matrix
	
	Vals[0][0]=tstep;	//Save the length of a timestep for plotting purposes	
	Vals[1][0]=tfinal;	//Save the final time for plotting purposes

	for(int i =1; i<n+1; i++)	//Loop through all of the bodies
	{
		Vals[0][i]=m[i-1];	//Save the mass of body[i]
		Vals[1][i]=l[i-1];	//Save the length of body[i]
	}
	
	//System Initialization
	horizontal_drop(inits,n);	//Set the initial conditions
	inits[n]=1;
	//Save the initial conditions in the solution matrix
	for(int r=0; r<2*n; r++)
	{
		X[0][r]=inits[r];
	}

	//Numerical Integration
	
	for(int t=1; t<tlen; t++)	//Loop through every timestep
	{
		
		RK_45(inits,tstep,n,m,l,II,Y,cut_off,Zs);	//Find the solution at that timestep
		for(int i = 0; i<2*n;i++)	//Loop through the solution
		{
			X[t][i]=Y[i];	//Save the solution at that timestep in the solution matrix
			inits[i]=Y[i];	//Use the solution as the initial conditions for the next timestep
		}
		if ( cudaSuccess != cudaGetLastError() )
   		 printf( "Error!\n" );
		
	}

	//Solution Output
	npy_save_double("mat.npy", fortran_order, 2, shape1, &X[0][0]);	//Output the solution
	npy_save_double("Vals.npy",fortran_order,2,shape2,&Vals[0][0]);	//Output values to find energy
	
	//Free memory
	free(inits);	//Initial conditions
	free(Y);	//Solution to each timestep
	free(m);
	free(l);
	free(II);
	free(Zs);

	return EXIT_SUCCESS;	//Program completed successfully
}
