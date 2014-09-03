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
void RecDCA(double Zs[], int n, int i, double AF[], int cut_off,double Xs[],int gpu, int data);


//Main function
int main()
{
	int data=1;
	int	n=0;
	std::ofstream timedata;
	std::ofstream timedata2;
	timedata2.open("cpuassembletime2.mtx");
	timedata.open("gpuassembletime2.mtx");


	//Variable Declarations
	//Variable Declarations
	//Initial conditions
  //List of joints between bodies NOTE: This joint list is not used in this version
	

	/////////////////////////////////////////////////////


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

	
	//Save the initial conditions in the solution matrix

	//myfile << "\n";
	float timeValue;
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
	
		double *Zs;
	double *m=(double*)malloc(sizeof(double)*n);
	double *Xs=(double*)malloc(sizeof(double)*n*5*5);
	//std::ofstream myfile;
	//std::ofstream myfile2;
	//myfile2.open("Vals.mtx");
  	//myfile.open ("output.mtx");
	//System Setup


	Zs = (double*)malloc(sizeof(double)*n*26*6);
		cudaEventRecord( beginEvent, 0 );
		
		RecDCA(Zs, n, 0,m,0,Xs,1,data);
	cudaEventRecord( endEvent, 0 );
	cudaEventSynchronize( endEvent );
	cudaEventElapsedTime( &timeValue, beginEvent, endEvent );
	timedata<< timeValue << "  ";
	
	cudaEventRecord(beginEvent,0);
		RecDCA(Zs,n,0,m,0,Xs,0,0);
	cudaEventRecord(endEvent,0);
	cudaEventSynchronize(endEvent);
	
	cudaEventElapsedTime( &timeValue, beginEvent, endEvent );
timedata2<<timeValue<<"  ";

 	if ( cudaSuccess != cudaGetLastError() )
    printf( "Error!\n" );
	std::cout<<"xxx"<<cudaGetLastError()<<"xxx"<<std::endl;
	std::cout << n << std::endl;
	free(Zs);
	free(m);
	free(Xs);
	}
	
	timedata2.close();
	timedata.close();
	return EXIT_SUCCESS;	//Program completed successfully
}	
