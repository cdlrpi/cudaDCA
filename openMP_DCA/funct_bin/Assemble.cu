//	Assemble.cu
//
//This file contains the function that assembles bodies

#include <iostream>
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
void Assemble(double Zs[], double Xs[],double nZs[], double nXs[], int len, int odd, int n)
{	
	//Variable Declarations
	//Both varaibles are temporary variables used for matrix operations
	double z1[6][6];
	double z2[6][6];
	double A[6][6];

	//Loop through every body in oldbds and newbds.  j is used to 
	//reference a body in oldbds and i to reference a body in newbds
#pragma omp parallel for
for(int i = 0, j=0; i<len-odd; i++, j+=2)
	{  
		for(int r=0; r<6; r++)	//Loop through every row
		{
			for(int c=0; c<6; c++)	//Loop through every column
			{
				z1[r][c]=Zs[c+r*n*26+i*52+26];	//Save b2.z11 into z1
				z2[r][c]=Zs[c+r*n*26+i*52+19];  //save b1.z22 into z2
			}
		}

		get_X(z1,z2,A);//Get the intermediate quantity X and save it in A
		invert_X(A);//Invert X and put it in A
		
		for(int r = 0; r<5; r++)	//Loop through every row
		{
			for(int c=0; c<5; c++)	//Loop through every column
			{
				nXs[c+r*5*len+i*5]=A[r][c];
					//Save Xinv in the new body corresponding to the 
			}									//two used to construct X
		}
			
		make_W(A,z1);	//Using X inverse, construct the intermediate quantity W and save it in A
		
		for(int r =0; r<6; r++)
		{
			for(int c=0; c<6; c++)
			{
				z1[r][c]=Zs[c+r*n*26+i*52+32]; //z212
				z2[r][c]=Zs[c+r*n*26+i*52+6];	//z112
			}
		}
		Mat66Mult(A,z1,z1);	//Perform A*b2.z12 and save the result in z1
		//z1 now holds W*b2.z12

		Mat66Mult(z2,z1,z1);	//Perform b1.z12*z1 and save the result in z1
		//z1 now holds b1.z12*W*b2.z12=z12
		
		for(int r = 0; r<6; r++)	//Loop through every row
		{
			for(int c=0; c<6; c++)	//Loop through every column
			{
				nZs[c+r*len*26+i*26+6]=z1[r][c];
				//Save the new z12 in the corresponding new body
				z1[r][c]= Zs[c+r*n*26+i*52+13];//z121
				z2[r][c]=Zs[c+r*n*26+i*52+6];//z112
		    }
		}
			
		Mat66Mult(A,z1,z1);	//Perform A*b1.z21 and store the result in z1
		//z1 now holds W*b1.z21

		Mat66Mult(z2,z1,z1);	//Perform b1.z12*z1 and store the result in z1
		//z1 now holds b1.z12*W*b1.z21
		
		for(int r = 0; r<6; r++)	//Loop through every row
		{
			for(int c=0; c<6; c++)	//Loop through every column
			{
				nZs[c+r*len*26+i*26]=Zs[c+r*26*n+i*52]-z1[r][c];	//Save the new z11 in newbds
				z1[r][c]=Zs[c+r*26*n+i*52+13];//z121
				z2[r][c]=Zs[c+r*26*n+i*52+39];//z221
			}
		}
		    	
		Mat66Mult(A,z1,z1);	//Perform A*b1.z21 and save the result in z1
		//z1 now holds W*b1.z21

		Mat66Mult(z2,z1,z1);	//Perform b2.z21*z1 and store the result in z1
		//z1 now holds b2.z21*W*b1.z21=z21
		
		for(int r = 0; r<6; r++)	//Loop through every row
		{
			for(int c=0; c<6; c++)	//Loop through every column
			{
				nZs[c+r*len*26+i*26+13]=z1[r][c];	//Save the new z21 in newbds
				z1[r][c]=Zs[c+r*26*n+i*52+32];//z212
				z2[r][c]=Zs[c+r*26*n+i*52+39];//z221
    		}
		}		
 
		Mat66Mult(A,z1,z1);	//Perform A*b2.z12 and store the result in z1
		//z1 now holds W*b2.z12

		Mat66Mult(z2,z1,z1);	//Perform b2.z21*z1 and store the result in z1
		//z1 now holds b2.z21*W*b2.z12
		
    	for(int r = 0; r<6; r++)	//Loop through every row
		{
			for(int c=0; c<6; c++)	//Loop through every column
			{
				nZs[c+r*len*26+i*26+19]=Zs[c+r*26*n+i*52+45]-z1[r][c];
    		}
		}
  		
		for(int r = 0; r<6; r++)	//Loop through every row
		{	
			z1[r][0]=Zs[r*26*n+i*52+25]-Zs[r*26*n+i*52+38];
			//Save b1.z23+b2.z13 into z1
			for(int c=0; c<6; c++)
			{
				z2[r][c]=Zs[c+r*26*n+i*52+6];
			}
		}
		
    	Mat61Mult(A,z1,A);	//Perform A*z1 and store the result in A
		//A now holds W*(b1.z23+b2.z13)=Y
    		
		Mat61Mult(z2,A,z1);	//Perform b1.z12*A and store the result in z1
		//z1 now holds b1.z12*Y
		
		for(int r = 0; r< 6; r++)	//Loop through every row
		{
			nZs[r*len*26+i*26+12]=Zs[r*n*26+i*52+12]-z1[r][0];
		
			//Save the new z13
			for(int c=0; c<6; c++)
			{
				z2[r][c]=Zs[c+r*26*n+52*i+39];
			}
		}
	
		Mat61Mult(z2,A,z1);	//Perform b2.z21*A and store the result in z1
		//z1 now holds b2.z21*Y
		
		for(int r=0; r< 6; r++)	//Loop through every row
		{
			nZs[r*len*26+i*26+25]=Zs[r*n*26+i*52+51]+z1[r][0];
				//Save the new z23 in newbds
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
				nZs[c+r*26*len+(len-1)*26]=Zs[c+r*26*n+(n-1)*26]; //z11
				nZs[c+r*26*len+(len-1)*26+6]=Zs[c+r*26*n+(n-1)*26+6]; //z12
				nZs[c+r*26*len+(len-1)*26+13]=Zs[c+r*26*n+(n-1)*26+13]; //z21
				nZs[c+r*26*len+(len-1)*26+19]=Zs[c+r*26*n+(n-1)*26+19]; //z22
			}
			nZs[r*26*len+(len-1)*26+12]=Zs[r*26*n+(n-1)*26+12];
			nZs[r*26*len+(len-1)*26+25]=Zs[r*26*n+(n-1)*26+25];
		}
	}
}
