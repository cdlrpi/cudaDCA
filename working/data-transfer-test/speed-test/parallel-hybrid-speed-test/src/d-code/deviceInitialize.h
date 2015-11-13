//	deviceInitialize.h
//
//This file contains the function used to find the initial zeta values of every body in the system

//Included Files
#include <stdio.h>
#include<cuda.h>

//Function Prototypes
//	Functions found in deviceFuncts.h
__device__ void MSsetup_gpu(double Minv[6][6],double S[6][6],int row, int col, double m, double I[],int Iindex);
__device__ void Fa_gpu(double S[6][6],double w, double r[], double m,int row, int col);
__device__ void zab_gpu(double r[],double z[6][6],double Minv[6][6],double S[6][6], int row, int col);
__device__ void zaa_gpu(double r[],double z[6][6],double Minv[6][6],double S[6][6], int row, int col);
__device__ void r02_gpu(double r[], double angle,double L);
__device__ void r01_gpu(double r[], double angle,double L);
__device__ void Mat61Mult_gpu(double A[6][6], double B[6][6], double C[6][6],int row, int col);
__device__ void Mat66Mult_gpu(double A[6][6], double B[6][6], double C[6][6],int row, int col);

//Initialize:
//	Function used to find the initial zeta values for every body in the system.
//		m is an array of the mass of each body in the system
//		l is an array of the length of every body in the system
//		I is an array holding the inertia tensors of all bodies in the system
//		Zetas is the array that will hold all of the zeta values for all bodies
//		n is the number of bodies
__global__ void Initialize_gpu(double state[],double m[], double l[],double I[],double Zetas[],int n)
{
    //Variable Declarations
	//	Variables to distinguish threads
	const int i1 = blockIdx.x;
	const int row = threadIdx.y;
	const int col = threadIdx.x;

	//	Indices to write solution
    const int Iindex = blockIdx.x*3+row*3*n+col;
    const int i11 = col+row*26*n+i1*26;
    const int i12=i11+6;
    const int i13=i12+6; //only col=0 can write
    const int i21=i13+1;
    const int i22=i21+6;
    const int i23=i22+6; //only col = 0 can write

	//	Temporary shared matrices for matrix operations
    __shared__ double z[6][6];
    __shared__ double Minv[6][6];
    __shared__ double S[6][6];
    __shared__ double r[3];
    __shared__ double q;
    __shared__ double w;

    int i = 0;	//counter
	q=0;
	w=0;

    while (i <= i1)	//Loop through the necessary amount of values
    {
        q+=state[i];	//Save the angle of the body
        w+=state[i+n];	//Save the angular velocity of the body
        i++;
    }
	__syncthreads();	//ensure w and q are found
   
    MSsetup_gpu(Minv,S,row,col,m[i1],I,Iindex);	//Set up the shifter and inverse mass matrix
    z[row][col]=0;	//reset z
    r01_gpu(r,q,l[i1]);	//Find the r01 position vector and save it in r
    zaa_gpu(r,z,Minv,S,row,col);	//Find z11 and store it in z
    Zetas[i11]=z[row][col];	//Save z11 into the output zeta array

    z[row][col]=0;	//reset z
    r02_gpu(r,q,l[i1]);	//Find the r02 position vector and save it in r
	zab_gpu(r,z,Minv,S,row,col);	//Find z12 and store it in z
    Zetas[i12]=z[row][col];	//Save z12 into the output zeta array

    r01_gpu(r,q,l[i1]);	//Find the r01 position vector and save it in r
    Fa_gpu(S,w,r,m[i1],row,col);	//Find the state dependent forces and save them in S
    Mat61Mult_gpu(Minv,S,z,row,col);	//Perform Minv*S and save the result in z

    if (col==0)	//Ensure only one column of threads enters the if
    {
        Zetas[i13]=z[row][col];	//Save z13 into the zeta matrix
    }

    z[row][col]=0;	//reset z
    MSsetup_gpu(Minv,S,row,col,m[i1],I,Iindex);	//Set up the shifter and inverse mass matrix
    r02_gpu(r,q,l[i1]);	//Find the r02 position vector
    zaa_gpu(r,z,Minv,S,row,col);	//Find z22 and store it in z
    Zetas[i22]=z[row][col];	//Save z22 in the zeta arary

    z[row][col]=0;	//reset z
    r01_gpu(r,q,l[i1]);	//Find the r01 position vector
    zab_gpu(r,z,Minv,S,row,col);	//Find z21 and store it in z
    Zetas[i21]=z[row][col];	//Save z21 in the zeta array
 
    r02_gpu(r,q,l[i1]);	//Find the r02 position vector 
    Fa_gpu(S,w,r,m[i1],row,col);	//Find the state dependent forces
    Mat61Mult_gpu(Minv,S,z,row,col);	//Perorm Minv*S and save the result in z

    if (col==0)	//Ensure only one column of threads enters the if
    {
        Zetas[i23]=z[row][col];	//save z23 in the zeta array
    }

	__syncthreads();	//Sync before finishing function
}
