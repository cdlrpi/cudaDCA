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
#include <string.h>
#include <math.h>
#include <cuda.h>
#include <fstream>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
//	created files
#include "d_code/deviceDisassemble.h"
#include "d_code/deviceAssemble.h"
#include "d_code/deviceInitialize.h"
#include "d_code/deviceFuncts.h"
#include "funct_bin/npy.h"
#include "funct_bin/classes.h"


//Function Prototypes
//	Function found in RK45.cu
void RK_45(double state[], double step, int n,double m[], double l[], double II[],double Y[],int cut_off,  std::vector<Body> body, std::vector<Forces> AF, double Zs[],float times[], int reps);

//	Functions found in Functs.cu
void pend_init(double m[], double l[], double II[],int n,double mass, double length);
void horizontal_drop(double x[],int n);
void set_up(double A[], double B[], double C[], int n , double h);

//	Functions found in Igaoibn;aeirn
void Initialize(double m[], double l[],double II[],double Zetas[],int n);
void printa(double A[], int n);
//Body::Body(int n);

// Functions found in findCutoff.cu
int findCutoff(int n, int accuracy);

using namespace std;	//namespace for output and input

//Main function
int main()
{
	int reps=40;
	int	n=0;
	int cut_off;
	float *times=(float*)malloc(sizeof(float)*reps);
	std::ofstream timedata;
	std::ofstream numbods;
	numbods.open("numbods_omp1-1.mtx");
	timedata.open("graph_ompdca1-1.mtx");
int numa =0;
	while(n<4000)
	{
		if(n<500)
		{
			n+=3;
		}
		else if( n<2000)
		{
			n+=10;
		}
		else if(n< 10000)
		{
			n+= 100;
		}
		else
		{
			n+=10000;
		}

//n+=5;
	cut_off = numa;
	

	std::cout<<"xxxx"<<cut_off<<"xxxx";
	//Variable Declarations
	//Variable Declarations
	double *inits;	//Initial conditions
	double *Y;	//Solution to each timestep
	double *m;
	double *l;
	double *II;  //List of joints between bodies NOTE: This joint list is not used in this version
	double *Zs;

	//std::ofstream myfile;
	//std::ofstream myfile2;
	//myfile2.open("Vals.mtx");
  	//myfile.open ("output.mtx");
	//System Setup
		vector<Body> bodies;
	vector<Forces> AF;
	int x=n;
	while(x != 1)
	{	
		Body b(x);
		Forces a(x);
		bodies.push_back(b);
		AF.push_back(a);
		if(x%2==0)
		{
			x/=2;
		}
		else
		{
			x++;
			x/=2;
		}
	}
	Body b(x);
	Forces a(x);
	bodies.push_back(b);
	AF.push_back(a);

	Zs = (double*)malloc(sizeof(double)*6*26*n);
	inits = (double*)malloc(sizeof(double)*2*n);	//Initial conditions are length 2*n
	Y = (double*)malloc(sizeof(double)*2*n);	//Timestep solution is length 2*n
	m = (double*)malloc(sizeof(double)*n);
	l = (double*)malloc(sizeof(double)*n);
	II = (double*)malloc(sizeof(double)*n*3*3);
	pend_init(m,l,II,n,1.0,1.0); //Initialize mass, length, and inertia of all bodies 

	Initialize(m, l, II, Zs, n);

	/////////////////////////////////////////////////////
	


	//Time Setup
	double tstep= 0.001; //Length of a timestep [s]
	double tfinal = 0.01; //Final time [s]
	int tlen = (int) floor(tfinal/tstep)+1;	//Number of timesteps

	//Matrix Output Setup
	//int shape1[2] = { tlen , 2*n }, fortran_order = 0;	//Shape of solution matrix
	//int shape2[2] = { 2 , n+1 };	//Shape of matrix holding information to calculate the energy
	//double Vals[2][n+1];	//Matrix holding information to calculate the energy
	

	typedef std::numeric_limits< double > dbl;

	
	std::cout.precision(dbl::digits10);
	//myfile2<<tstep<<"  ";
	//Vals[0][0]=tstep;	//Save the length of a timestep for plotting purposes	
	//Vals[1][0]=tfinal;	//Save the final time for plotting purposes


	//System Initialization
	horizontal_drop(inits,n);	//Set the initial conditions
	

		RK_45(inits,tstep,n,m,l,II,Y,numa,bodies, AF,Zs,times,reps);	//Find the solution at that timestep



	//	for(int i = 0; i<2*n;i++)	//Loop through the solution
		//{
			
		//	inits[i]=Y[i];	//Use the solution as the initial conditions for the next timestep
			//myfile << inits[i]<<"  ";
	//	}
		//myfile << "\n";
		
//	}
	float timeValue=0;
	for(int i =0; i<reps; i++)
	{
		timeValue+=times[i];
	}
	timeValue/=reps;
	timedata<< timeValue << "  ";
	numbods<<n<<"  ";
 	if ( cudaSuccess != cudaGetLastError() )
    printf( "Error!\n" );
	std::cout << n << std::endl;
	//Solution Output
	//npy_save_double("Vals.npy",fortran_order,2,shape2,&Vals[0][0]);	//Output values to find energy
	
	//Free memory

	
	//myfile.close();
	//myfile2.close();
	//std::cout<<n<<std::endl;
	}
	
	numbods.close();
	timedata.close();
	return EXIT_SUCCESS;	//Program completed successfully
}	
