//	cudaDCA.cu
//
//This file contains the recursive DCA function, and the function that is used to invoke DCA and
//interperate the results.

//Included Files

#include <iostream>
#include <vector>
#include "classes.h"
//Function Prototypes
//	Functions found in this file
void RecDCA(std::vector<Body> bodies, std::vector<Forces> AF, int n,int cut_off);

//	Functions found in Init_setup.cu
void CudaInitialize(double m[], double l[], double I[], double x[], int n, double Zs[]);

//	Functions found in Assemble_setup.cu
void cudaAssemble(double Zs[],double Xs[], int num, double nZs[], double nXs[], int odd, int newlen);

//	Functions found in Disassemble_setup.cu
void cudaDisassemble(double OldAF[], double Zs[], double Xs[],double nZs[], double nXs[], int odd, int morelen, int lesslen, double AF[]);

//	Functions found in Assemble.cu
void Assemble(double Zs[], double Xs[],double nZs[], double nXs[], int len, int odd, int n);

//	Functions found in Disassemble.cu
void Disassemble(double lessZs[], double lessXs[],double moreZs[], double moreXs[], double oldAs[] ,double newAs[], int num, int odd);

//	Functions found in SolveBCs.cu
void solve_BCs(double Zs[], double Xs[], double AF[]);

void printa(double A[], int n);
void printm(double A[6][6]);

void Update_Properties(double bodyZetas[],double nZetas[], int n, double state[], double m[], double l[], double II[]);
//DCAhelp:
//	Function that prepares the list of bodies for DCA and finds the final state vector
//		state is the state of the system at that timestep
//		bs is a list of bodies used for initialization
//		js is a list of joints
//		n is the number of bodies
//		Y is the array where the final velocities and accelerations are stored
void DCAhelp(double state[], double m[], double l[], double I[],int n, double Y[],int cut_off, std::vector<Body> bodies , std::vector<Forces> AF, double *Zs, float times[], int reps)
{
	//Create the list that will hold all acceleration and force values for all bodies
	

	Update_Properties(Zs,bodies[0].Zs,n,state,m,l,I);

	//CudaInitialize(m,l,I, state, n, Zs);	//Initialize the bodies, finding all zeta values
	
	//Pass the list of bodies to DCA and return the accelerations 
	//and forces of both handles of every body in the list


	
	cudaEvent_t beginEvent;
	cudaEvent_t endEvent;
	cudaEventCreate( &beginEvent );
	cudaEventCreate( &endEvent );
	for(int i= 0; i<reps; i++)
	{
		cudaEventRecord( beginEvent, 0 );
	RecDCA(bodies,AF, n, cut_off);
	cudaEventRecord( endEvent, 0 );
	cudaEventSynchronize( endEvent );
	cudaEventElapsedTime( &times[i], beginEvent, endEvent );
	}
	Y[n]=AF[0].af[8*n];	//For a pendulum, the fist acceleration value is in A[2][0]

	for(int i = n+1, j=2; i<n*2; i++, j+=2)	//Loop through the acceleration matrix
	{
		Y[i]=AF[0].af[2*4*n+2*j]-AF[0].af[2*4*n+2*(j-1)]; //Find and save all generalized accelerations
	}

	for(int i = 0; i<n; i++)	//Loop through the state vector
	{
		Y[i]=state[i+n];	//Save the velocities
	}

	//Free memory
		
}

//RecDCA:
//	Function used to solve for the velocty and acceleration of the list of bodies at
//	the current timestep.  This is a recursive function that continues to call itself
//	until there is a single body left.  Once at this point the accelerations and forces
//	are found using the boundary conditions of a pendulum.  These values are then returned
//	to the previous level of recursion which then finds the new accelerations and forces
//	for the disassembled bodies.  This continues until all bodies are disassembled, ultimately
//	returning the forces and accelerations at both handles of every body in the system.  These
//	results are intererated by DCAhelp (above) to obtain the actual generalized accelerations.  
//		bodies is the list of bodies
//		n is the number of bodies
//		i is the level of recursion
//		AF is the array in which the accelerations and forces at the handles of the bodies 
//			will be stored.
void RecDCA(std::vector<Body> bodies, std::vector<Forces> AF, int n,int cut_off)
{
	int *nums = (int*)malloc(sizeof(int)*bodies.size());
	int x = n;
	int newlen;
	int odd;
	int i =0;
	while(x!=1)
	{
		odd=0;
		nums[i]=x;
		if(x%2==0)
		{
			newlen = (int) (x/2);
		}
		else
		{
			odd =1;
			newlen=(int)((x+1)/2);
		}
		/*
		if(i<cut_off)
		{
			cudaAssemble(bodies[i].Zs,bodies[i].Xs, x, bodies[i+1].Zs, bodies[i+1].Xs , odd, newlen);	//Assemble the bodies, storing them in newbds			
		}
		else
		{
 	*/
			Assemble(bodies[i].Zs,bodies[i].Xs,bodies[i+1].Zs,bodies[i+1].Xs, newlen,odd, x);	//Assemble the bodies, storing them in newbds
		//}
		i++;
		x=newlen;
	}
	nums[i]=x;

	//Call the DCA function again to return the accelerations and forces of the new bodies
	
	solve_BCs(bodies[i].Zs,bodies[i].Xs, AF[i].af);	//Solve the boundary conditions and find the acceleratins and forces
		//Knowing the accelerations and forces of the new bodies, the new bodies can be disassembled
		//again, finding the accelerations and forces of the old bodies.
	while(i!=0)
	{
		newlen = nums[i];
		x=nums[i-1];
		if(x%2==0)
		{
			odd = 0;
		}
		else
		{
			odd = 1;
		}
	
		//if(i-1<cut_off)
		//{
			
		//	cudaDisassemble(AF[i].af, bodies[i-1].Zs,bodies[i-1].Xs , bodies[i].Zs,bodies[i].Xs, odd, x, newlen, AF[i-1].af);
			
		//}
	
		//else

	//	{
			Disassemble(bodies[i].Zs,bodies[i].Xs,bodies[i-1].Zs,bodies[i-1].Xs,AF[i].af, AF[i-1].af, newlen,odd);
		//}

		i--;
	}
		 
	
	
}
		
			

