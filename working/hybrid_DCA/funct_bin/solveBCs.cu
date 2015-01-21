//	solveBCs.cu
//
//This file contains the function used solve for the boundary conditions of a pendulum

//Included Files
#include "classes.h"
#include <iostream>

//Function Prototypes
//	Functions found in Functs.cu
void matmult61(double A[6][6], double B[6], double C[6]);

//solve_BCs:
//	Function used to solve for the acceleration and force at each handle of a pendulum.
//	This function does assume that the first handle is pinned and the second is free,
//	as is the case for this simulation.
//		bodies is the single body that is a product of all assembled bodies
//		AF is where the accelerations and forces will be stored
void solve_BCs(Body *bodies, double AF[])
{
	//Variable declarations
	double Fc1[6];	//Constraint force on handle 1
	double Fc2[6];	//Constraint force on handle 2
	double A1[6];	//Acceleration of handle 1
	double A2[6];	//Acceleration of handle 2
	double temp[6];	//Temporary matrix used for matrix operations
	double val;		//Temporary value
	double M[3][4];	//Matrix used to solve a system of linear equations

	//This loop fills the M matrix with the correct values in z11 and z13 in order to solve
	//for the force at handle 1.  This can be done because handle 1 of a pendulum is pinned,
	//therefor no translational acceleration is allowed in the joint, and no rotational forces
	//are allowed in the joint
	for(int c = 0; c<3; c++)	//Loop through 3 columns
	{
		for(int r = 0; r<3; r++)	//Loop through 3 rows
		{
			M[r][c] = bodies[0].z11[r+3][c+3];	//Save the correct values in z11

			if (c==0)	//If the column is 0
			{
				M[r][3] = -1*bodies[0].z13[r+3];	//Save the correct values in z13
			}
		}
	}
	
	//This loop solves the system of linear equations created in the loop above
	//by using gaussian elimination to turn the leftmost 3x3 matrix in M into the 
	//identity matrix
	for(int s=0; s<3;s++)	//Loop through 3 rows
	{
		val = M[s][s];	//Set the temporary value so it is not overwritten durring the loop
		for(int r=0; r<4;r++)	//Loop through every column
		{	
			M[s][r]=M[s][r]/val;	//Divide every element in the row by the first non zero element
		}
		for(int j =0; j<3; j++)	//Loop through 3 rows
		{	
			val=M[j][s];	//Set the temporary value so it is not overwritten durring the loop
			
			if( s!=j)	//If the current iteration is not on the diagonal
			{
				for(int l=0; l<4; l++)	//Loop through every column
				{	
					//Subtract the rows to produce zeros below the leading 1
					M[j][l]=M[j][l]-(M[s][l]*val);	
				}
			}
		}
	}
	//After completing the loop above, the last column of M holds the translational 
	//constraint force on the first handle of the body

	for(int s=0; s<6;s++)	//Loop through every row
	{
		//Set the constraint force on the second handle to 0 becuase this end is free
		Fc2[s] = 0;

		if (s>=3)	//If the row is greater than 3
		{
			Fc1[s] = M[s-3][3];	//Save the translational force from M into Fc1
		}
		else
		{
			Fc1[s] = 0;	//Set the rest of Fc1 to 0 because there is no rotational force
		}
	}
	
	matmult61(bodies[0].z11, Fc1, temp);	//Perform z11*Fc1 and save the result in temp
	
	for(int s = 0; s<6; s++)	//Loop through every row
	{
		A1[s]=temp[s]+bodies[0].z13[s];	//Find and save A1
	}
	
	matmult61(bodies[0].z21, Fc1, temp);	//Perform z21*Fc1 and save the result in temp

	for(int s = 0; s<6; s++)	//Loop through every row
	{
		A2[s]=temp[s]+bodies[0].z23[s];	//Find and save A2
	}
	
	for(int r =0; r<6; r++)	//Loop through every row
	{
		AF[r*4]=A1[r];	//Save A1 into AF
		AF[r*4+1]=Fc1[r];	//Save Fc1 into AF
		AF[r*4+2]=A2[r];	//Save A2 into AF
		AF[r*4+3]=Fc2[r];	//Save Fc2 into Af
	}
	
}

//The function below is used for debugging and printing purposes only and can be
//removed at any time
void printM(double M[3][4])
{
	for(int i = 0; i < 3; i++)
	{
		for( int j = 0; j<4;j++)
		{
			std::cout<<M[i][j]<<"\t";
		}
		std::cout<<std::endl;
	}
}
