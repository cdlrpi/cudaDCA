//
//	3D DCA Driver
//This file invokes all of the necessary function calls to prepare
//and simulate a compound pendulum system through the use of the
//recursive DCA algorithm.	The majority of this algorithm is run
//on the gpu.  Output is created in a format that is
//readable in python for answer checking and graphing purposes.

//Included Files
#include <malloc.h>
#include <iostream>
#include <math.h>
#include <math.h>
#include <fstream>
#include <limits>

#define z11 0
#define z12 6
#define z13 12
#define z21 18
#define z22 24
#define z23 30

#define nz11 0
#define nz12 6
#define nz13 12
#define nz21 13
#define nz22 19
#define nz23 25
//Function Prototypes
//	Function found in RK45.cu
void Initialize(double Mass[],double Inertia[], double Zetas[], int n, bool Active_DOF[], double Body_Placement[],double Body_Vectors[],int dof_index[]);
void init_2d_pend(bool Active_DOF[],double Body_Vectors[], double Mass[], double Inertia[],int n, double Body_Placement[], double speeds[]);
void printZeta(double[],int,int,int,int);
void update(double Mass[], double Inertia[], double Init_Zetas[], double Zetas[], int n, double Body_Vectors[], double Speeds[], double DCMs[]);
void kinematics(bool Active_DOF[], double Coords[], double DCMs[], double speeds[], double omegas[], int dof_index[],int n);
void RK_45(double step, int n,bool Active_DOF[],double Coords[], double Speeds[], double initZetas[],double Mass[], double Inertia[],int DOF, double Y[], double Ydot[], int dof_index[], int cut_off, double Body_Vectors[]);
//Main function
int main()
{
	int cut_off=0;
	int	n=2;
	int DOF=0;
	std::ofstream output;
	output.open("output.mtx");


	//Variable Declarations
	double *Body_Placement = new double[n*6];	//Initial conditions
	double *Body_Speeds = new double[n*6];
	double *Mass = new double[n];
	double *Inertia = new double[3*n];
	bool *Active_DOF = new bool[n*6];		//Solution to each timestep
	double *Body_Vectors = new double[n*6];
	int *dof_index = new int[n];

	init_2d_pend(Active_DOF,Body_Vectors,Mass,Inertia,n,Body_Placement,Body_Speeds);

	for(int i =0; i<6*n; i++)
	{
		if(Active_DOF[i])
		{
			DOF++;
		}
	}
	
	double *Coords = new double[DOF];
	double *Speeds = new double[DOF];
	double *initZetas = new double[n*36*6];

	
	for(int i =0, j=0; i<n*6; i++)
	{
		if(Active_DOF[i])
		{
			Coords[j]=Body_Placement[i];
			Speeds[j]=Body_Speeds[i];
			j++;
		}
	}

	Initialize(Mass, Inertia, initZetas, n, Active_DOF, Body_Placement,Body_Vectors,dof_index);
	//printZeta(initZetas,n,0,36,z11);
	//double *DCMs = new double[n*3*3];
	//double *omegas= new double[n*3];
	//kinematics(Active_DOF, Coords, DCMs, Speeds, omegas, dof_index, n);
	//update(Mass, Inertia, initZetas, Zetas, n, Body_Vectors, omegas, DCMs);
		
	//Time Setup
	double tstep= 0.001; //Length of a timestep [s]
	double tfinal = 0.005; //Final time [s]
	int tlen = (int) floor(tfinal/tstep)+1;	//Number of timesteps

	double *Y = new double[DOF];
	double *Ydot = new double[DOF];
	for(int t=1; t<tlen; t++)	//Loop through every timestep
	{
		RK_45(tstep,n,Active_DOF,Coords,Speeds,initZetas,Mass, Inertia,DOF,Y,Ydot,dof_index,cut_off,Body_Vectors);	//Find the solution at that timestep

		for(int i = 0; i<DOF;i++)	//Loop through the solution
		{
			Coords[i]=Y[i];
			Speeds[i]=Ydot[i];
		}
	}

	delete[] Body_Placement; 
	delete[] Body_Speeds;
	delete[] Mass;
	delete[] Inertia; 
	delete[] Active_DOF;
	delete[] Body_Vectors;
	delete[] dof_index;
	delete[] Y;
	delete[] Ydot;
	delete[] Coords;
	delete[] Speeds;
	delete[] initZetas;

	return EXIT_SUCCESS;	//Program completed successfully
}



	
void init_2d_pend(bool Active_DOF[],double Body_Vectors[], double Mass[], double Inertia[],int n, double Body_Placement[], double speeds[])
{
	double m =1;
	//double r = 0.05;
	double l =1;
/*
	double I1 = .5*m*r*r;
	double I2 = (m/12)*((3*r*r)+(l*l));
	double I3 = I2;
*/
	double I1 = (m/12);
	double I2 = (m/12);
	double I3 = (m/12); 
	 
	for(int i =0; i<6*n; i++)
	{
		Body_Placement[i]=0;
		speeds[i]=0;
		Body_Vectors[i]=0;
	}
	for(int i =0; i<n; i++)
	{

		Active_DOF[i*6+2]=1;

		Mass[i]=1.0;

		Body_Vectors[i*6]=-l/2;
		Body_Vectors[i*6+3]=l/2;
		
		Inertia[i*3]=I1;
		Inertia[i*3+1]=I2;
		Inertia[i*3+2]=I3;
	}
	//Body_Placement[2]=0.3;
		
}		
void printZeta(double Zetas[],int n, int body,int len,int zeta)
{
	std::cout<<"\n\n\n";
	for(int r =0; r<6; r++)
	{
		for(int c =0; c<6; c++)
		{
			std::cout<<Zetas[r*n*len+body*len+zeta+c]<<'\t';	
		}
		std::cout<<std::endl;
	}
}		
