//	Initialize.cu
//
//This file contains the function used to find the initial zeta values of every body in the system

//Included Files
#include "classes.h"
#include <stdio.h>
#include <iostream>

//Function Prototypes
//	Functions found in DCAfuncts.cu
void MSsetup(float Minv[6][6],float S[6][6],float m, float I[3][3]);
void Fa(float S[6][6],float w, float r[], float m);
void r02(float r[], float angle,float L);
void r01(float r[], float angle,float L);
void Mat66Mult(float A[6][6], float B[6][6], float C[6][6]);
void Mat61Mult(float A[6][6], float B[6][6], float C[6][6]);
void zaa(float r[],float z[6][6],float Minv[6][6],float S[6][6]);
void zab(float r[],float z[6][6],float Minv[6][6],float S[6][6]);


//Initialize:
//	Function used to find the initial zeta values for every body in the system.
//		state is an array of the state of the system at that timestep
//		oldbds is the list of bodies that do not yet have zeta values
//		newbds is the list of bodies where the new zeta values will be saved
//		n is the number of bodies
void Initialize(float state[],InitBody *oldbds,Body *newbds,int n)
{
    //Variable Declarations
    float z[6][6];
    float Minv[6][6];
    float S[6][6];
    float rr[3];
    float q;
    float w;
	
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
