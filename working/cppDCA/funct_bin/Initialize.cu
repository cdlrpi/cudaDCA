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
void r02(double r[], double angle,double L);
void r01(double r[], double angle,double L);
void Mat66Mult(double A[6][6], double B[6][6], double C[6][6]);
void Mat61Mult(double A[6][6], double B[6][6], double C[6][6]);
void zaa(double r[],double z[6][6],double Minv[6][6],double S[6][6]);
void zab(double r[],double z[6][6],double Minv[6][6],double S[6][6]);


//Initialize:
//	Function used to find the initial zeta values for every body in the system.
//		state is an array of the state of the system at that timestep
//		oldbds is the list of bodies that do not yet have zeta values
//		newbds is the list of bodies where the new zeta values will be saved
//		n is the number of bodies
void Initialize(double state[],InitBody *oldbds,Body *newbds,int n)
{
    //Variable Declarations
    double z[6][6];
    double Minv[6][6];
    double S[6][6];
    double rr[3];
    double q;
    double w;
	
	for(int i =0; i< n; i++)	//Loop through every body
	{	
		q=0;	//Initialize q to 0
		w=0;	//Initialize w to 0
		int j =0;	//counter

       	while (j <= i)	//Loop through the necessary amount of values
    	{
	        q+=state[j];	//Save the angle of the body
	        w+=state[j+n];	//Save the angular velocity of the body
	        j++;
    	}
		
       	MSsetup(Minv,S,oldbds[i].m,oldbds[i].I);	//Set up the shifter and inverse mass matrix
    
		for(int r =0; r< 6; r++)	//Loop through every row
		{
			for(int c=0; c<6; c++)	//Loop through every column
			{
				z[r][c]=0;	//Initialize z to 0
			}
		}
	
		r01(rr,q,oldbds[i].l);	//Find the r01 position vector and save it in rr
		zaa(rr,z,Minv,S);	//use r01, Minv, and S to find z11
		
		for(int r =0; r< 6; r++)	//Loop through every row
		{
			for(int c=0; c<6; c++)	//Loop through every column
			{
				newbds[i].z11[r][c]=z[r][c];	//Save z11
				z[r][c]=0;	//reset z
			}
		}	

		r02(rr,q,oldbds[i].l);	//Find the r02 position vector and save it in rr
		zab(rr,z,Minv,S);	//Use r02, Minv, and S to find z12

		for(int r =0; r< 6; r++)	//Loop through every row
		{
			for(int c=0; c<6; c++)	//Loop throughe every column
			{
				newbds[i].z12[r][c]=z[r][c];	//Save z12
			}
		}

		r01(rr,q,oldbds[i].l);	//Find the r01 position vector and save it in rr
		Fa(S,w,rr,oldbds[i].m);	//Find the state dependent force and save it in S
		Mat61Mult(Minv,S,z);	//Perform Minv*S and save the result in z
		//z now holds z13

		for(int r =0; r< 6; r++)	//Loop throughe every row
		{
	   		newbds[i].z13[r]=z[r][0];	//Save z13
		}
    
		for(int r =0; r< 6; r++)	//Loop through every row
		{
			for(int c=0; c<6; c++)	//Loop through every column
			{	
				z[r][c]=0;	//Reset z
			}
		}

	    MSsetup(Minv,S,oldbds[i].m,oldbds[i].I); //Set up the shifter and inverse mass matrix
	    r02(rr,q,oldbds[i].l);	//Find the r02 position vector and save it in rr
	    zaa(rr,z,Minv,S);	//use r02, Minv, and S to find z22

		for(int r =0; r< 6; r++)	//Loop through every row
		{
			for(int c=0; c<6; c++)	//Loop through every column
			{
				newbds[i].z22[r][c]=z[r][c];	//Save z22
				z[r][c]=0;
			}
		}
   
		r01(rr,q,oldbds[i].l);	//Find the r01 position vector and save it in rr
		zab(rr,z,Minv,S);	//Use r01, Minv, and S to find z21

		for(int r =0; r< 6; r++)	//Loop through every row
		{
			for(int c=0; c<6; c++)	//Loop through every column
			{
				newbds[i].z21[r][c]=z[r][c];	//Save z21
			}
		}

		r02(rr,q,oldbds[i].l);	//Find the r02 position vector and save it in rr
		Fa(S,w,rr,oldbds[i].m);	//Find the state dependent force and save it in S
		Mat61Mult(Minv,S,z);	//Perform Minv*S and store the result in z

	    for(int r =0; r< 6; r++)	//Loop through every row
	    {
			newbds[i].z23[r]=z[r][0];	//Save z23
	    }
	}	//End loop through bodies
}
