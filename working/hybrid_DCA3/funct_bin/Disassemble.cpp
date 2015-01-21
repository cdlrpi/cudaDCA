//	Disassemble.cu
//
//This file contains the function that disassembles bodies.
#define z111 0
#define z112 6
#define z113 12
#define z121 13
#define z122 19
#define z123 25
#define z211 26
#define z212 32
#define z213 38
#define z221 39
#define z222 45
#define z223 51
//Included Files

#include <iostream>

//Function Prototypes
//	Functions found in DCAfuncts.cu
void Mat61Mult(double A[6][6], double B[6][6], double C[6][6]);
void make_W(double Xinv[6][6], double D[6][6]); 

//Disassemble:
//	Function used to disassemble bodies and solve for the forces and accelerations at the joints.
//		lessbds is a list of the assembled bodies
//		morebds is a list of the disassembled bodies
//		oldAs is a list of the accelerations and forces for the assembled bodies
//		newAs is a list of the accelerations and forces for the disassembled bodies
//		num is the number of assembled bodies
//		odd is equal to 1 if there is an odd number of disassembed bodies and equal to 
//			0 if there is an even number.	
void Disassemble(double lessZs[], double lessXs[],double moreZs[], double moreXs[], double oldAs[] ,double newAs[], int num, int odd)
{	
	//Variable Declarations
	//All variables are temporary variables used to perform matrix operations	
	double Vals[6][6];
	double temp1[6][6];
	double temp2[6][6];
	double temp3[6][6];
	double A[6][6];
	int morelen=(int)((num*2)-odd);

	for(int i=0; i<num-odd; i++)	//Loop through every assembled body
	{	

		for(int r=0; r<6; r++)	//Loop through every row
		{
			for(int c=0; c< 6; c++)	//Loop through every column
			{
				temp1[r][c]=0;	//initialize temp1 to 0
				temp2[r][c]=0;	//initialize temp2 to 0
				
				if(c<4)	//If the column is less than 4
				{
					Vals[r][c]=oldAs[c+r*(num)*4+i*4];	//Save the forces and accelerations for the
				}										//current body in Vals
    							
				if(c==0) //If the column is 0
				{
					temp1[r][c]=moreZs[r*morelen*26+i*52+z123]-moreZs[r*morelen*26+i*52+z213]; //Save b1.z23+b2.z13 in temp1
				}
			}
		}
		
		
		for(int r=0; r<6; r++)	//Loop through every row
		{
    		temp2[r][0]=Vals[r][3];	//Save b2.Fc2 into column 0 of temp2
			for(int c=0; c<6; c++)
			{
				temp3[r][c]=moreZs[c+r*morelen*26+i*52+z212];
			}
   		}

    	Mat61Mult(temp3,temp2,temp2);	//perform the operation b2.z12*temp2 and save
		//temp 2 now holds b2.z12*b2.Fc2			//the result in temp2
			
		for(int r=0; r<6; r++)	//Loop through every row
		{
			temp2[r][0]=-1*temp2[r][0];	//make temp2 negative
			//temp2 now holds -1*b2.z12*b2.Fc2

			temp2[r][0]=temp2[r][0]+temp1[r][0];		//Add temp2 and temp1 and store the result 
			//temp2 now holds -1*z212*F2c2+z123-z213	//in temp2

			A[r][0]=Vals[r][1];	//Save b1.Fc1 in the first row of A
			for(int c=0; c<6; c++)
			{
				temp3[r][c]=moreZs[c+r*morelen*26+i*52+z121];
			}
    	}

		Mat61Mult(temp3,A,A);	//perform b1.z21*A and store the resul in A
		//A now holds b1.z21*b1.Fc1
		
		for(int r=0; r<6; r++)	//Loop through every row
		{
			for(int c=0; c<6; c++)	//Loop through every column
			{
				if(c==0)	//If the column is 0
				{
    				temp2[r][c]=A[r][c]+temp2[r][c];	//add A to temp2
					//temp2 now holds b1.21*b1.c1-b2.12*b2.c2+b1.23-b2.13
				}
    	
    			if(r<5 && c<5)	//If the row and colum is less than 5
    			{
    				A[r][c]=lessXs[c+r*num*5+i*5];	//Load X inverse into A
    			}
			}		
		}
		
    	make_W(A,temp1);	//Use X inverse to create the intermediate matrix W
    	//A now holds W

    	Mat61Mult(A,temp2,A);	//Perform the operation A*temp2 and store the result in A
		//A now holds W*(b1.21*b1.c1-b2.12*b2.c2+b1.23-b2.13)

		for(int r=0; r<6; r++)	//Loop through every row
		{
    		A[r][0]=-1*A[r][0];	//A now holds -1*(W*((z121*F1c1)-(z212*F2c2)+z123-z213))
			//A is now equal to b1.Fc2

    		newAs[r*(num*2-odd)*4+i*8+3]=A[r][0];	//Save b1.Fc2 into the newAs list
    		temp2[r][0]=Vals[r][1];	//Load b1.Fc1 into temp2
    	}
    	
    	Mat61Mult(temp3,temp2,temp1);	//Perform the operation b1.z21*temp2 and store the
		//temp1 now holds b1.z21*b1.Fc1			//result in temp1
    	for(int r = 0; r<6; r++)
		{
			for(int c=0; c<6; c++)
			{
				temp3[r][c]=moreZs[c+r*morelen*26+i*52+z122];
			}
		}
		Mat61Mult(temp3,A,temp2);	//Perform the operation b1.z22*A and store the 
		//temp2 now holds b1.z22*b1.Fc2		//result in temp2

		for(int r=0; r<6; r++)	//Loop through every row
		{
    		temp1[r][0]=temp1[r][0]+temp2[r][0];	//Add temp2 to temp1
        	temp2[r][0]=moreZs[r*morelen*26+i*52+z123];	//save b1.z23 into temp2
    		temp2[r][0]=temp1[r][0]+temp2[r][0];	//temp2 now holds b1.A2
        	newAs[r*(num*2-odd)*4+i*8+2]=temp2[r][0];	//Save b1.A2 in the newAs list
    		A[r][0]=A[r][0]*-1;	//A now holds b2.Fc1
        	newAs[r*(num*2-odd)*4+i*8+5]=A[r][0];	//Save b2.Fc1 in the newAs list
			temp2[r][0]=Vals[r][3];	//Load b2.Fc2 into temp2
			for(int c=0; c<6; c++)
			{
				temp3[r][c]=moreZs[c+r*morelen*26+i*52+z211];
			}
   		}
   
    	Mat61Mult(temp3,A,A);	//Perform the operation b2.z11*A and store the
		//A now holds b2.z11*b2.Fc1			//result in A
		for(int r=0; r<6; r++)
		{
			for(int c=0; c<6; c++)
			{
				temp3[r][c]=moreZs[c+r*morelen*26+i*52+z212];
			}
		}
    	Mat61Mult(temp3,temp2,temp2);	//Perform the operation b2.z12*temp2 and store
		//temp2 now holds b2.z12*b2.Fc2				//the result in temp2
		
		for(int r=0; r<6; r++)	//Loop through every row
		{
			temp1[r][0]=A[r][0]+temp2[r][0]+moreZs[r*morelen*26+i*52+z213];	//temp1 now holds b2.A1    	
        	newAs[r*(num*2-odd)*4+i*8+4]=temp1[r][0];	//Save b2.A1 in the newAs list
			
			//Save all the values in oldAs into the correct spot in newAs
        	newAs[r*(num*2-odd)*4+i*8]=Vals[r][0];
        	newAs[r*(num*2-odd)*4+i*8+1]=Vals[r][1];
        	newAs[r*(num*2-odd)*4+i*8+6]=Vals[r][2];
        	newAs[r*(num*2-odd)*4+i*8+7]=Vals[r][3];
    	}

	}	//End loop through bodies
	//std::cout<<"a"<<std::endl;
	//If there is an odd number of bodies, the forces and accelerations from
	//the last body in lessbds must be added to the newAs list
	if(odd==1)
	{
		for(int r=0; r<6; r++)	//Loop through every row
		{
			//Save the forces and accelerations from the last body into the new list
			newAs[(r+1)*(num*2-odd)*4-4]=oldAs[(r+1)*num*4-4];
			newAs[(r+1)*(num*2-odd)*4-3]=oldAs[(r+1)*num*4-3];
			newAs[(r+1)*(num*2-odd)*4-2]=oldAs[(r+1)*num*4-2];
			newAs[(r+1)*(num*2-odd)*4-1]=oldAs[(r+1)*num*4-1];
		}
	}
}
