//	Assemble.cu
//
//This file contains the function that assembles bodies

//Included Files
#include "classes.h"

//Function Prototypes
//	Functions found in DCAfuncts.cu
void Mat66Mult(double A[6][6], double B[6][6], double C[6][6]);
void Mat61Mult(double A[6][6], double B[6][6], double C[6][6]);
void get_X(double z1[6][6], double z2[6][6], double D[6][6]);
void invert_X(double X[6][6]);
void make_W(double Xinv[6][6], double D[6][6]); 
void printm(double a[6][6]);
//Assemble:
//	Function used to assemble a list of bodies into a list of bodies that is 
//	half the size of the original list. To accomplish this, the list of old bodies
//	is cycled through looking at two bodies at a time.  These two bodies are assembled
//	into one new body.  The counter then moves on to look at the next two old bodies
//  If the original list has an odd number of bodies, the last body is ignored during
//	assembly and is added on to the new list of bodies at the end of the function.
//		oldbds is the old list of bodies
//		newbds is the new list of assembled bodies
//		len is the length of oldbds
//		odd is 1 if oldbds has an odd number of bodies and 0 if it has an even number
void Assemble(Body *oldbds,Body *newbds, int len, int odd)
{	
	//Variable Declarations
	//Both varaibles are temporary variables used for matrix operations
	double z1[6][6];
	double A[6][6];

	//Loop through every body in oldbds and newbds.  j is used to 
	//reference a body in oldbds and i to reference a body in newbds
	for(int i = 0, j=0; i<len-odd; i++, j+=2)
	{  
		for(int r=0; r<6; r++)	//Loop through every row
		{
			for(int c=0; c<6; c++)	//Loop through every column
			{
				z1[r][c]=oldbds[j+1].z11[r][c];	//Save b2.z11 into z1
			}
		}

		get_X(z1,oldbds[j].z22,A);//Get the intermediate quantity X and save it in A
		invert_X(A);//Invert X and put it in A
		
		for(int r = 0; r<5; r++)	//Loop through every row
		{
			for(int c=0; c<5; c++)	//Loop through every column
			{
				newbds[i].Xinv[r][c]=A[r][c];	//Save Xinv in the new body corresponding to the 
			}									//two used to construct X
		}
			
		make_W(A,z1);	//Using X inverse, construct the intermediate quantity W and save it in A
		
		Mat66Mult(A,oldbds[j+1].z12,z1);	//Perform A*b2.z12 and save the result in z1
		//z1 now holds W*b2.z12

		Mat66Mult(oldbds[j].z12,z1,z1);	//Perform b1.z12*z1 and save the result in z1
		//z1 now holds b1.z12*W*b2.z12=z12
		
		for(int r = 0; r<6; r++)	//Loop through every row
		{
			for(int c=0; c<6; c++)	//Loop through every column
			{
				newbds[i].z12[r][c]=z1[r][c];	//Save the new z12 in the corresponding new body
		    }
		}
			
		Mat66Mult(A,oldbds[j].z21,z1);	//Perform A*b1.z21 and store the result in z1
		//z1 now holds W*b1.z21

		Mat66Mult(oldbds[j].z12,z1,z1);	//Perform b1.z12*z1 and store the result in z1
		//z1 now holds b1.z12*W*b1.z21
		
		for(int r = 0; r<6; r++)	//Loop through every row
		{
			for(int c=0; c<6; c++)	//Loop through every column
			{
	  			newbds[i].z11[r][c]=oldbds[j].z11[r][c]-z1[r][c];	//Save the new z11 in newbds
			}
		}
		    	
		Mat66Mult(A,oldbds[j].z21,z1);	//Perform A*b1.z21 and save the result in z1
		//z1 now holds W*b1.z21

		Mat66Mult(oldbds[j+1].z21,z1,z1);	//Perform b2.z21*z1 and store the result in z1
		//z1 now holds b2.z21*W*b1.z21=z21

		for(int r = 0; r<6; r++)	//Loop through every row
		{
			for(int c=0; c<6; c++)	//Loop through every column
			{
   	  			newbds[i].z21[r][c]=z1[r][c];	//Save the new z21 in newbds
    		}
		}		
 
		Mat66Mult(A,oldbds[j+1].z12,z1);	//Perform A*b2.z12 and store the result in z1
		//z1 now holds W*b2.z12

		Mat66Mult(oldbds[j+1].z21,z1,z1);	//Perform b2.z21*z1 and store the result in z1
		//z1 now holds b2.z21*W*b2.z12
	
    	for(int r = 0; r<6; r++)	//Loop through every row
		{
			for(int c=0; c<6; c++)	//Loop through every column
			{
    			newbds[i].z22[r][c]=oldbds[j+1].z22[r][c]-z1[r][c];	//save the new z22 into newbds
    		}
		}
  
		for(int r = 0; r<6; r++)	//Loop through every row
		{
			z1[r][0]=oldbds[j].z23[r]-oldbds[j+1].z13[r];	//Save b1.z23+b2.z13 into z1
		}

    	Mat61Mult(A,z1,A);	//Perform A*z1 and store the result in A
		//A now holds W*(b1.z23+b2.z13)=Y
    		
		Mat61Mult(oldbds[j].z12,A,z1);	//Perform b1.z12*A and store the result in z1
		//z1 now holds b1.z12*Y

		for(int r = 0; r< 6; r++)	//Loop through every row
		{
			newbds[i].z13[r]= oldbds[j].z13[r]-z1[r][0];	//Save the new z13
		}
	
		Mat61Mult(oldbds[j+1].z21,A,z1);	//Perform b2.z21*A and store the result in z1
		//z1 now holds b2.z21*Y

		for(int r=0; r< 6; r++)	//Loop through every row
		{
			newbds[i].z23[r]= oldbds[j+1].z23[r]+z1[r][0];	//Save the new z23 in newbds
		}
	}	//End loop through  bodies
	
	//If there is and odd number of oldbds, the list can not be cut directly in half.
	//Because of this, the last body in oldbds must be added to the end of newbds.
	if(odd ==1)
	{	
		for(int r=0; r<6; r++)	//Loop through every row
		{
			for(int c=0; c<6; c++)	//Loop through every column
			{
				newbds[len-1].z11[r][c]=oldbds[2*len-2].z11[r][c];	//Save z11
				newbds[len-1].z12[r][c]=oldbds[2*len-2].z12[r][c];	//Save z12
				newbds[len-1].z21[r][c]=oldbds[2*len-2].z21[r][c];	//Save z21
				newbds[len-1].z22[r][c]=oldbds[2*len-2].z22[r][c];	//Save z22
			}
			newbds[len-1].z13[r]=oldbds[2*len-2].z13[r];	//Save z13
			newbds[len-1].z23[r]=oldbds[2*len-2].z23[r];	//Save z23
		}
	}
}
