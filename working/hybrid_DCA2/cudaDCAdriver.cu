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
#include "funct_bin/classes.h"
#include "d_code/deviceDisassemble.h"
#include "d_code/deviceAssemble.h"
#include "d_code/deviceInitialize.h"
#include "d_code/deviceFuncts.h"
#include "funct_bin/npy.h"
#include <math.h>
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
	int	n=0;
	int cut_off;
	std::ofstream timedata;
	std::ofstream numbods;
	numbods.open("numbods.mtx");
	timedata.open("graph2_init_on_cpu.mtx");
	for(int numa = 4 ; numa<12; numa++)
	{
	n=0;
	while(n<80000)
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
		else
		{
			n+=10000;
		}
	int x = n;
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
	//std::cout<<x<<std::endl;
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
	double *II;  //List of joints between bodies NOTE: This joint list is not used in this version

	//std::ofstream myfile;
	//std::ofstream myfile2;
	//myfile2.open("Vals.mtx");
  	//myfile.open ("output.mtx");
	//System Setup
	
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

	//Matrix Output Setup
	//int shape1[2] = { tlen , 2*n }, fortran_order = 0;	//Shape of solution matrix
	//int shape2[2] = { 2 , n+1 };	//Shape of matrix holding information to calculate the energy
	//double Vals[2][n+1];	//Matrix holding information to calculate the energy
	

	typedef std::numeric_limits< double > dbl;

	
	std::cout.precision(dbl::digits10);
	//myfile2<<tstep<<"  ";
	//Vals[0][0]=tstep;	//Save the length of a timestep for plotting purposes	
	//Vals[1][0]=tfinal;	//Save the final time for plotting purposes
	cudaEvent_t beginEvent;
	cudaEvent_t endEvent;
	cudaEventCreate( &beginEvent );
	cudaEventCreate( &endEvent );

	//System Initialization
	horizontal_drop(inits,n);	//Set the initial conditions
	
	//Save the initial conditions in the solution matrix

	//myfile << "\n";
	cudaEventRecord( beginEvent, 0 );
	//Numerical Integration
	for(int t=1; t<tlen; t++)	//Loop through every timestep
	{
		RK_45(inits,tstep,n,m,l,II,Y,cut_off);	//Find the solution at that timestep
		for(int i = 0; i<2*n;i++)	//Loop through the solution
		{
			
			inits[i]=Y[i];	//Use the solution as the initial conditions for the next timestep
			//myfile << inits[i]<<"  ";
		}
		//myfile << "\n";
		
	}
	cudaEventRecord( endEvent, 0 );
	cudaEventSynchronize( endEvent );
	float timeValue;
	cudaEventElapsedTime( &timeValue, beginEvent, endEvent );
	timedata<< timeValue << "  ";
	numbods<<n<<"  ";
 	if ( cudaSuccess != cudaGetLastError() )
    printf( "Error!\n" );
	std::cout << n << std::endl;
	//Solution Output
	//npy_save_double("Vals.npy",fortran_order,2,shape2,&Vals[0][0]);	//Output values to find energy
	
	//Free memory
free(inits);
	free(II);
	free(m);
	free(l);
	free(Y);
	
	//myfile.close();
	//myfile2.close();
	//std::cout<<n<<std::endl;
	}
	timedata<<"\n";
	numbods<<"\n";
	}
	
	numbods.close();
	timedata.close();
	return EXIT_SUCCESS;	//Program completed successfully
}	
