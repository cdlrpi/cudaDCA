//	DCA.cu
//
//This file contains the recursive DCA function, and the function that is used to invoke DCA and
//interperate the results

//Included Files
#include "classes.h"
#include <iostream>

//Function Prototypes
//	Functions found in this file
void RecDCA(Body *bodies, int n, int i, float AF[]);

//	Functions found in Initialize.cu
void Initialize(float state[],InitBody *oldbds,Body *newbds,int n);

//	Functions found in Assemble.cu
void Assemble(Body *oldbds,Body *newbds, int len, int odd);

//	Functions found in Disassemble.cu
void Disassemble(Body *lessbds,Body *morebds, float oldAs[] ,float newAs[], int num, int odd);

//	Functions found in SolveBCs.cu
void solve_BCs(Body *bodies, float AF[]);

//DCAhelp:
//	Function that prepares the list of bodies for DCA and finds the final state vector
//		state is the state of the system at that timestep
//		bs is a list of bodies used for initialization
//		js is a list of joints
//		n is the number of bodies
//		Y is the array where the final velocities and accelerations are stored
void DCAhelp(float state[], InitBody *bs, Joint *js,int n, float Y[])
{
	//Create the list that will hold all acceleration and force values for all bodies
	float *AF = (float*) malloc(sizeof(float)*n*4*6);
	float A[6][2*n];	//Create the matrix where only the accelerations will be stored
	Body *bodies = new Body[n];	//Create the list of bodies that will be use in DCA

	Initialize(state,bs, bodies, n);	//Initialize the bodies, finding all zeta values

	//Pass the list of bodies to DCA and return the accelerations 
	//and forces of both handles of every body in the list
	RecDCA(bodies, n, 0, AF);	

	for(int r = 0; r<6; r++)	//Loop through every row
	{
		for(int c = 0; c<2*n;c++)	//Loop through every column
		{
			A[r][c]=AF[c*2+r*4*n];	//Save only the accelerations from AF into A
		}
	}

	Y[n]=A[2][0];	//For a pendulum, the fist acceleration value is in A[2][0]

	for(int i = n+1, j=2; i<n*2; i++, j+=2)	//Loop through the acceleration matrix
	{
		Y[i]= A[2][j]-A[2][j-1];	//Find and save all generalized accelerations
	}

	for(int i = 0; i<n; i++)	//Loop through the state vector
	{
		Y[i]=state[i+n];	//Save the velocities
	}

	//Free memory
	free(AF);
	delete[] bodies;		
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
void RecDCA(Body *bodies, int n, int i, float AF[])
{
	if (n==1)	//If there is only 1 body
	{
		solve_BCs(bodies, AF);	//Solve the boundary conditions and find the acceleratins and forces
	}
	else	//If there is more than 1 body
	{	
		int newlen;	//New number of bodies after assembly
		int odd = 0;	//Flag to keep track of the parity of the length of the list of bodies
		Body *newbds;	//Create the list of newly assembled bodies

		if(n % 2 == 0)	//If there is an even number of bodies
		{
			newlen = (int) (n/2);	//The new number of bodies will be half the original number
		}
		else	//If there is an odd number of bodies
		{
			newlen = (int)((n+1)/2);	//The new number of bodies will be half the original number
										//rounded down, plus 1
			odd = 1;	//odd is set to 1 because there are an odd number of bodies
		}
		
		newbds = new Body[newlen];	//Initialize the new list to the new length
		Assemble(bodies, newbds, newlen,odd);	//Assemble the bodies, storing them in newbds
		
		//Create a list of accelerations and forces of the new bodies.
		float *AFo=(float*)malloc(sizeof(float)*6*newlen*4);
		
		//Call the DCA function again to return the accelerations and forces of the new bodies
		RecDCA(newbds,newlen,i+1 ,AFo);	
		
		//Knowing the accelerations and forces of the new bodies, the new bodies can be disassembled
		//again, finding the accelerations and forces of the old bodies.
		Disassemble(newbds,bodies,AFo, AF, newlen,odd);
		
 		//Free memory
		free(AFo);
		delete[] newbds;
	}
}
		
			

