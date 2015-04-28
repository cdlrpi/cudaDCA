//	Initialize.cu
//
//This file contains the function used to find the initial zeta values of every body in the system

//Included Files
#include "classes.h"
#include <stdio.h>
#include <iostream>

//Function Prototypes
//	Functions found in DCAfuncts.cu
void MSsetup(double Minv[6][6],double S[6][6],double m, double I[3][3]);
void Fa(double S[6][6],double w, double r[], double m);
void r02(double r[], double L);
void r01(double r[], double L);
void Mat66Mult(double A[6][6], double B[6][6], double C[6][6]);
void Mat61Mult(double A[6][6], double B[6][6], double C[6][6]);
void zaa(double r[],double z[6][6],double Minv[6][6],double S[6][6]);
void zab(double r[],double z[6][6],double Minv[6][6],double S[6][6]);
void printm(double A[6][6]);
void printa(double A[], int n);

//Initialize:
//	Function used to find the initial zeta values for every body in the system. n2 is up.
//		state is an array of the state of the system at that timestep
//		oldbds is the list of bodies that do not yet have zeta values
//		newbds is the list of bodies where the new zeta values will be saved
//		n is the number of bodies
void Initialize(double m[], double l[],double II[],double Zetas[],int n)
{
    //Variable Declarations
	//	Variables to distinguish threads


	//	Temporary shared matrices for matrix operations
    double z[6][6];
    double Minv[6][6];
    double S[6][6];
    double rr[3];
	double I[3][3];
  
	
	for(int j=0; j<n; j++)
	{
	
   		for(int r =0; r<3; r++)
		{
			for(int c=0; c<3; c++)
			{
				I[r][c]=II[j*3+r*3*n+c];
			}
		}
				
    	MSsetup(Minv,S,m[j],I);	//Set up the shifter and inverse mass matrix
		
		for(int r =0; r<6; r++)
		{
			for(int c=0; c<6; c++)
			{
    			z[r][c]=0;	//reset z					
			}
		}

    	r01(rr,l[j]);	//Find the r01 position vector and save it in r
		
    	zaa(rr,z,Minv,S);	//Find z11 and store it in z
		
		for(int r =0; r<6; r++)
		{
			for(int c=0; c<6; c++)
			{
    			Zetas[c+r*26*n+j*26]=z[r][c];	//Save z11 into the output zeta array
				z[r][c]=0;	//reset z
			}
		}
 
    	r02(rr,l[j]);	//Find the r02 position vector and save it in r
		zab(rr,z,Minv,S);	//Find z12 and store it in z
		for(int r =0; r<6; r++)
		{
			for(int c=0; c<6; c++)
			{
    			Zetas[c+r*26*n+j*26+6]=z[r][c];	//Save z12 into the output zeta array
				z[r][c]=0;
			}
			Zetas[12+r*26*n+j*26]=0;
		}
    
		MSsetup(Minv,S,m[j],I);	//Set up the shifter and inverse mass matrix
		r02(rr,l[j]);	//Find the r02 position vector
		zaa(rr,z,Minv,S);	//Find z22 and store it in z
		for(int r =0; r<6; r++)
		{
			for(int c=0; c<6; c++)
			{

    			Zetas[c+r*26*n+j*26+19]=z[r][c];	//Save z22 in the zeta arary
				z[r][c]=0;	//reset z
			}
		}
   
		r01(rr,l[j]);	//Find the r01 position vector
		zab(rr,z,Minv,S);	//Find z21 and store it in z

		for(int r =0; r<6; r++)
		{
			for(int c=0; c<6; c++)
			{
    			Zetas[c+r*26*n+j*26+13]=z[r][c];	//Save z21 in the zeta array
			}
			Zetas[r*26*n+j*26+25]=0;
		}

	}
}
