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

	//User Input and Display
	cout<<"\n";	//clear screen

	//If there are no arguments, default to outputting an mtx file
	cout<<"Number of bodies:  ";
	cin>>n;
	flag = mtx_file;

	//Ask user for final time
	cout<<"\n\n"<<"Final time:\t";
	cin>>tfinal;
	tlen = (int) floor(tfinal/tstep)+1;	//length of time matrix

	//Begin Initializing system and determining number of gpu assemblies
	cout<<"\n\nDetermining most efficient configuration...\n\n";
	numa=findCutoff(n,20);	//Find number of gpu assemblies

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
	
	//Ask user to press a key to begin
	cout<<numa<<" assemblies will be performed on the gpu\n\nPress enter to begin";
	c=getchar();
	c=getchar();
	cout<<"\nNumber of bodies:\t"<<n<<"\nFinal time:\t"<<tfinal<<"\nNumber of assemblies performed on gpu:\t"<<numa<<"\n\n\t\t";

	ofstream myfile;
	ofstream myfile2;
	myfile2.open("Body_Values.mtx");
  	myfile.open ("Solution.mtx");
	
	//Set maximum output precision
	std::cout.precision(16);
	myfile2<<tstep<<"  ";

	//Matrix Output Setup
	for(int i =1; i<n+1; i++)	//Loop through all of the bodies
	{
		myfile2<<m[i-1]<<"  ";	//Record the mass of every body
	}
	myfile2<<"\n"<<tfinal<<"  ";//Record the final time
	for(int i =1; i<n+1; i++)	//Loop through all of the bodies
	{
		myfile2<<l[i-1]<<"  ";	//Record the length of every body
	}

	//Save the initial conditions in the solution matrix
	for(int r=0; r<2*n; r++)
	{
		myfile << inits[r]<< "  ";
	}
	myfile << "\n";

	for(int t=1; t<tlen; t++)	//Loop through every timestep
	{
		RK_45(inits,tstep,n,m,l,II,Y,numa,Zs);	//Find the solution at that timestep
		for(int i = 0; i<2*n;i++)	//Loop through the solution
		{
			inits[i]=Y[i];	//Use the solution as the initial conditions for the next timestep
			myfile << inits[i]<<"  ";
		}
		if((((float)t)/((float)tlen))>=percheck)	//Check for progress
		{
			cout<<percent<<"%\n\t\t";	//Print progress
			percheck+=.1;
			percent+=10;
		}
	myfile << "\n";
	}

	cout<<"100%\n\t\t";

	//Check for errors in integration and complete output
	if ( cudaSuccess != cudaGetLastError() )
	{
	 	cout<<"\nAn error occurred during integration, No output generated\n";
		return 0;
	}
	else
	{	
		myfile.close();
		myfile2.close();
		cout<<"\nIntegration successful! Type \"python plot.py\" to plot results\n";
	}
	free(inits);	//Initial conditions
	free(Y);	//Solution to each timestep
	free(m);
	free(l);
	free(II);
	free(Zs);

	return EXIT_SUCCESS;	//Program completed successfully
}
