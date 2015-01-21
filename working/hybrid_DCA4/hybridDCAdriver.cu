//
//	Hybrid DCA Driver
//
//This file invokes all of the necessary function calls to prepare
//and simulate a compound pendulum system through the use of the
//recursive DCA algorithm.	Much of this algorithm is run
//on the gpu.  Output can be chosen using the flags -p -m or -b
//
//-p generates a more accurate output file that can only be read in python.
//-m generates a file with printed numbers that can be read by any software
//-b generates both
//
//When creating a python file, there may be memory issues for large numbers of bodies.
//Because of this, the number of bodies for a python file is limited to 100

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
#include "d_code/deviceDisassemble.h"
#include "d_code/deviceAssemble.h"
#include "d_code/deviceInitialize.h"
#include "d_code/deviceFuncts.h"
#include "funct_bin/npy.h"

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
	cout<<string( 100, '\n' );	//clear screen
	if(argc==2)	//argument check
	{
		if(argv[1][1]=='p' || argv[1][1]=='b')
		{
			//Ask user for number of bodies
			cout<<"Number of bodies (less than 100):\t";
			cin>>n;	
			
			//Check number of bodies
			if(n>100)
			{
				cout<<"/n/nIncorrect number of bodies!";
				return 0;
			}

			//Set the flag to keep track of output
			if(argv[1][1]=='p')
			{
				flag = python_file;
			}
			else
			{
				flag = both;
			}
		}
		else if(argv[1][1]=='m')
		{
			flag = mtx_file;	//Set flag
			
			//Ask user for number of bodies
			cout<<"Number of bodies:  ";
			cin>>n;
		}
	}
	//If there are too many arguments, leave the program
	else if(argc>2)
	{
		cout<<"Too many arguments!\nProgram Terminated!";
		return 0;
	}
	//If there are no arguments, default to outputting an mtx file
	else if(argc==1)
	{
		cout<<"Number of bodies:  ";
		cin>>n;
		flag = mtx_file;
	}

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
	cout<<string(100,'\n')<<"Number of bodies:\t"<<n<<"\nFinal time:\t"<<tfinal<<"\nNumber of assemblies performed on gpu:\t"<<numa<<"\n\n\t\t";
	
	//Python output file
	if(flag == python_file)
	{

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

		//Save the initial conditions in the solution matrix
		for(int r=0; r<2*n; r++)
		{
			X[0][r]=inits[r];
		}

		//Numerical Integration
		for(int t=1; t<tlen; t++)	//Loop through every timestep
		{
		
			RK_45(inits,tstep,n,m,l,II,Y,numa,Zs);	//Find the solution at that timestep
			for(int i = 0; i<2*n;i++)	//Loop through the solution
			{
				X[t][i]=Y[i];	//Save the solution at that timestep in the solution matrix
				inits[i]=Y[i];	//Use the solution as the initial conditions for the next timestep
			}
			
			if((((float)t)/((float)tlen))>=percheck)	//Check if enough progress was made
			{
				cout<<percent<<"%\n\t\t";	//Print progress
				percheck+=.1;
				percent+=10;
			}
		
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
			//Solution Output
			npy_save_double("pySol.npy", fortran_order, 2, shape1, &X[0][0]);	//Output the solution
			npy_save_double("bodyVals.npy",fortran_order,2,shape2,&Vals[0][0]);	//Output values to find energy
			
			cout<<"\nIntegration successful! Type \"ipython plotpy.py\" to see results\n";
		}
	}

	//Mtx file output
	else if(flag == mtx_file)
	{
		ofstream myfile;
		ofstream myfile2;
		myfile2.open("bodyVals.mtx");
	  	myfile.open ("mtxSol.mtx");
		
		//Set maximum output precision
		typedef numeric_limits< double > dbl;
		std::cout.precision(15);
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
			cout<<"\nIntegration successful! Type \"ipython plotmtx.py\" to see results\n";
		}
	}

	//Both files output
	else if(flag == both)
	{

		//Mtx file setup
		std::ofstream myfile;
		std::ofstream myfile2;
		myfile2.open("bodyVals.mtx");
	  	myfile.open ("mtxSol.mtx");
		
		//Set maximum output precision
		typedef std::numeric_limits< double > dbl;
		std::cout.precision(15);
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

		//Python file output
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

		//Save the initial conditions in the solution matrix
		for(int r=0; r<2*n; r++)
		{
			X[0][r]=inits[r];
		}
	
		for(int t=1; t<tlen; t++)	//Loop through every timestep
		{
			RK_45(inits,tstep,n,m,l,II,Y,numa,Zs);	//Find the solution at that timestep

			for(int i = 0; i<2*n;i++)	//Loop through the solution
			{
				inits[i]=Y[i];	//Use the solution as the initial conditions for the next timestep
				myfile << inits[i]<<"  ";
				X[t][i]=Y[i];	//Save the solution at that timestep in the solution matrix
			}
			if((((float)t)/((float)tlen))>=percheck)	//Check progress
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
			npy_save_double("pySol.npy", fortran_order, 2, shape1, &X[0][0]);	//Output the solution
			npy_save_double("bodyVals.npy",fortran_order,2,shape2,&Vals[0][0]);	//Output values to find energy
			cout<<"\nIntegration successful! Type \"ipython plotmtx.py\" to see results of mtx file, or \"ipython plotpy.py\" to see results of python file\n";
		}
	}

	//Free memory
	free(inits);	//Initial conditions
	free(Y);	//Solution to each timestep
	free(m);
	free(l);
	free(II);
	free(Zs);

	return EXIT_SUCCESS;	//Program completed successfully
}
