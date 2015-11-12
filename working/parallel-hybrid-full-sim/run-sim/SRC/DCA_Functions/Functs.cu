//	Functs.cu
//
//This file contains functions used to perform operations that are not specific to the DCA algorithm

//Included Files

#include <iostream>

//matmult61:
//	This function is used to multiply a 6x6 matrix and a 6x1 matrix.
//		A is the 6x6 matrix
//		B is the 6x1 matrix
//		C is the matrix in which the solution is saved
void matmult61(double A[6][6], double B[6], double C[6])
{	
	double j[6];	//Temporary matrix

	for(int r = 0; r<6; r++)	//Loop through every row
	{
		j[r]=0;	//Initialize the elements of j to 0
		
		//Loop through all of the columns, taking the dot product of the 
		//row of A and the column of B and storing it in j
		for(int c = 0; c<6; c++)
		{
			j[r]+= (A[r][c]*B[c]);
		}
	}

	//Store the solution in C
	for(int r=0; r<6; r++)
	{
		C[r]=j[r];
	}
}


//pend_init:
//	This function is used to initialize the bodies of the system to the given mass and length.
//	All bodies are given the same properties
//		bs is a list of the bodies
//		n is the number of bodies
//		mass is the mass of the bodies
//		length is the length of the bodies	
void pend_init(double m[], double l[], double I[],int n,double mass, double length)
{
	int k=0;	//counter
	double II;	//Inertia

	while(k<n)	//Loop through all of the bodies
	{
		m[k]=mass;	//set the mass
		l[k]=length;	//set the length
		II= mass*(pow(length,2)/12.0);	//Find the inertia

		//Set the inertia tensor
		for(int r =0; r<3; r++)
		{
			for(int c = 0; c<3; c++)
			{
				if(r==c)
				{
					I[c+3*n*r+k*3]= II;
				}
			}
		}
		k++;
	}
}

//horizontal_drop:
//	Function used to set the initial conditions to that of dropping the pendulum from a completely
//	horizontal position.  Note that the angle of the first link in the pendulum is measured with
//	respect to the vertical.  All other angles are measured with respect to the previous link
//		X is the array that the initial conditions will be saved in
//		n is the number of bodies
void horizontal_drop(double x[],int n)
{
	int k = 0;	//counter
	x[k]=3.14159/2;	//Set the angle of the first link to a horizontal position
	k++;
	
	//Loop through the rest of the initial conditions, setting all other angles
	//and velocities to 0
	while(k<2*n)
	{
		x[k]=0;
		k++;
	}
}

//arycpy:
//	Function used to copy array B into array A
//		A is the array where the copy will be stored
//		B is the array to be copied
//		n is the number of values in the array
void arycpy(double A[],double B[],int n)
{
	int i=0;	//counter

	//Loop through the arrays saving every element in B into A
	while (i<n)
	{
		A[i]=B[i];
		i++;
	}
}

//arycpy2:
//	Function used to copy array B into the second half of array A.
//	Array A is assumed to be twice as long as B
//		A is the array where the copy will be stored
//		B is the array to be copied
//		n is the number of values in the array
void arycpy2(double A[], double B[], int n)
{
	int i = 0;	//counter

	//Loop through the arrays saving every element in B to the correct space in A
	while(i<n)
	{
		A[i+n]=B[i];
		i++;
	}
}

//arycpy3:
//	Function used to copy the second half of array B into array A.
//	Array B is assumed to be twice as long as A
//		A is the array where the copy will be stored
//		B is the array to be copied
//		n is the number of values in the array
//Function to copy the second half of array B into array A
void arycpy3(double A[], double B[], int n)
{
	int i = 0;	//counter

	//Loop through the arrays saving every element in B to the correct space in A
	while(i<n)
	{
		A[i]=B[i+n];
		i++;
	}
}

//set_up:
//	Function used to set up the qs and qdots in the numerical integrator
//		A is q or qdot
//		B is qa or qdota where a is between 1 and 4
//		C is where the solution is stored
//		n is the number of bodies
//		h is half of the length of one timestep
void set_up(double A[], double B[], double C[], int n , double h)
{
	int i = 0;	//counter

	//Loop through the elements in the arrays saving the solution in C
	while (i<n)
	{
		C[i]=A[i]+(h*B[i]);
		i++;
	}
}

//get_final_q:
//	Function used to get the final q to be used in the next timestep.
//		x is q
//		h is have of the length of one timestep
//		a is qdot1
//		b is qdot2
//		c is qdot3
//		d is qdot4
//		R is where the solution is stored
//		n is the number of bodies
void get_final_q(double x[], double h, double a[], double b[], double c[], double d[], double R[], int n)
{
	int i = 0;	//counter

	//Loop through the arrays solving for the position at the next timestep and saving it in R
	while(i<n)
	{
		R[i]= x[i] + ((h/3.0)*(a[i]+(2*b[i])+(2*c[i])+d[i]));
		i++;
	}
}

//get_final_qdot:
//	Function used to get the final qdot to be used in the next timestep.
//		x is qdot
//		h is have of the length of one timestep
//		a is qddot1
//		b is qddot2
//		c is qddot3
//		d is qddot4
//		R is where the solution is stored
//		n is the number of bodies
void get_final_qdot(double x[], double h, double a[], double b[], double c[], double d[], double R[], int n)
{
	int i = 0;	//counter

	//Loop through the arrays solving for the position at the next timestep and saving it in R
	while(i<n)
	{
		R[i+n]= x[i] + ((h/3.0)*(a[i]+(2*b[i])+(2*c[i])+d[i]));
		i++;
	}
}

/*
//Joint::Joint()
//	Joint initialization function
//	The joint is assumed to be that of a 2d pendulum
//	NOTE: This function is not used in this version
Joint::Joint()
{
	//int P[]={0,0,1,0,0,0};
	//int D[6][5]={{1,0,0,0,0},{0,1,0,0,0},{0,0,0,0,0},{0,0,1,0,0},{0,0,0,1,0},{0,0,0,0,1}};
}
*/

//////////////////////////////////////////////////////////////////////////////////////////////////////
//Functions below this line are used for debugging purposes only and can be removed at any time
void printm(double A[6][6])
{	
	for(int r = 0;r<6;r++)
	{
		for(int c = 0; c < 6; c++)
		{

			std::cout<<A[r][c]<<"\t";
		}
		std::cout<<std::endl;
	}
	std::cout<<std::endl<<std::endl;
}
void printx(double A[5][5])
{	
	for(int r = 0;r<5;r++)
	{
		for(int c = 0; c < 5; c++)
		{

			std::cout<<A[r][c]<<"\t";
		}
		std::cout<<std::endl;
	}
	std::cout<<std::endl<<std::endl;
}
void printa(double A[], int x)
{
	for(int i = 0; i<x;i++)
	{
		std::cout<<A[i]<<"  ";
	}
	std::cout<<std::endl<<std::endl;
}
