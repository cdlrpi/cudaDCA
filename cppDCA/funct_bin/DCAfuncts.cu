//	DCAfuncts.cu
//
//This file contains all of the functions that are used during initialization, 
//assembly, and disassembly.

//Function Prototypes
//	All function called in this file are written in this file
void transpose(float S[6][6]);
void Mat66Mult(float A[6][6], float B[6][6], float C[6][6]);
void Mat61Mult(float A[6][6], float B[6][6], float C[6][6]);
void makeS(float S[6][6],float r[]);
void set_D(float D[6][6]);
void set_Dt(float D[6][6]);

//get_X:
//	Function used to find the intermediate quantity X for a body
//	Considering the two bodies used to create X as b1 and b2:
//		z1 is a matrix containing z11 for b2
//		z2 is z22 for b1
//		D is the matrix in which X will be stored
void get_X(float z1[6][6], float z2 [6][6], float D[6][6])
{		
	for(int r = 0; r<6; r++)	//Loop through all 6 rows
	{
		for(int c=0; c<6; c++)	//Loop through all 6 columns
		{
    		z1[r][c]=z1[r][c]+z2[r][c];	//Add z1 and z2 together and store the solution in z1
		}
	}
	//Completing the above loop, z1 now holds b1.z22+b2.z11

    set_D(D);	//put the D matrix corresponding to a joint in this simulation in D
    Mat66Mult(z1,D,z1);	//Perform z1*D and store the solution in z1
	//z1 now holds (b1.z22+b2.z11)*D

    set_Dt(D);	//put the transpose of the D matrix in D
    Mat66Mult(D,z1,D); //Perform D*z1 and store the solution in D
	//D now holds Dt*(b1.z22+b2.z11)*D=X	
}

//set_D:
//	Function used to set the matrix that defines the subspace of prohibited motion.
//		D is the matrix in which the prohibited motion matrix will be stored.  
void set_D(float D[6][6])
{
	for(int r = 0; r<6; r++)	//Loop through all 6 rows
	{
		for(int c=0; c<6; c++)	//Loop through all 6 columns
		{
    		D[r][c]= 0;	//Set the D matrix to 0
		}
	}    
	
	//Set the D matrix according to the 2D joints used in this simulation
    D[0][0]=1;
    D[1][1]=1;
    D[3][2]=1;
    D[4][3]=1;
    D[5][4]=1;
}

//set_Dt:
//	Function used to set the transpose of the matrix that defines the subspace of prohibited motion.
//		D is the matrix in which the transposed prohibited motion matrix will be stored.  
void set_Dt(float D[6][6])
{
    for(int r = 0; r<6; r++)	//Loop through all 6 rows
	{
		for(int c=0; c<6; c++)	//Loop through all 6 columns
		{
    		D[r][c]= 0;	//Set the D matrix to 0;
		}
	}
   
	//Set the D matrix according to the 2D joints used in this simulation
    D[0][0]=1.0;
    D[1][1]=1.0;
    D[2][3]=1.0;
    D[3][4]=1.0;
    D[4][5]=1.0;
}

//invert_X:
//	Function used to invert the X matrix.  This function takes advantage of the nature of this 
//	2D simulation by only inverting the 2x2 matrix inside of X that will effect the solution.
//	This was done for efficiency of the simulation as well as ease of implementation in this case.
//		X is the intermediate quantity X; X inverse will be stored in X, replacing it
void invert_X(float X[6][6])
{
	//Variable Declarations
	float det;	//The determinant of the 2x2 matrix
	float Xinv[2][2];	//The inverse of the 2x2 matrix

	//Matrix inversion
	det=(X[2][2]*X[3][3])-(X[2][3]*X[3][2]);	//Find the determinant of the 2x2 matrix

	//Create the inverse matrix based on the trick for 2x2 matrices 
	Xinv[0][0]=X[3][3]/det;
	Xinv[0][1]=-1*X[2][3]/det;
	Xinv[1][0]=-1*X[3][2]/det;
	Xinv[1][1]=X[2][2]/det;

	//Save the Solution
	for(int c = 0; c<2; c++)	//Loop through 2 columns
	{
		for(int r=0; r<2;r++)	//Loop through 2 rows
		{
			X[r+2][c+2]=Xinv[r][c];	//Save X inverse in X
		}
	}
}

//make_W:
//	Function used to create the intermediate quantity W
//		Xinv is Xinverse;  W will be stored here
//		D is an empty matrix used for matrix calculations
void make_W(float Xinv[6][6], float D[6][6])                    
{   
	for(int c =0; c< 5; c++)	//Loop through 5 columns
	{
		for(int r=0; r<5; r++)	//Loop through 5 rows
		{
			//Set the values in row and column 6 to 0
    		if(r==5 || c==5)	
			{
     			Xinv[r][c]=0;	
			}
		}
	}

	set_Dt(D);	//Put the D transpose matrix in D
    Mat66Mult(Xinv,D,Xinv); //Perform Xinv*D and store the solution in Xinv
	//Xinv now holds Xinv*Dt
	
    set_D(D);	//Put the D matrix in D
    Mat66Mult(D,Xinv,Xinv);	//Perform D*Xinv and store the solution in Xinv
	//Xinv now holds D*Xinv*Dt=W
}

//zaa:
//	Function used to calculate z11 and z22 during initialization
//		r is the necessary position vector (r01 for z11 and r02 for z22)
//		z is the matrix in which the solution will be stored
//		Minv is the inverse mass matrix
//		S is an empty matrix that is made into the shifter matrix corresponding to the r vector
void zaa(float r[],float z[6][6],float Minv[6][6],float S[6][6])
{
    makeS(S,r);	//Create the shifter matrix
   	Mat66Mult(Minv,S,z);	//Perform Minv*S and store the solution in z
	//z now holds Minv*S

    transpose(S);	//Transpose the shifter matrix
    Mat66Mult(S,z,z);	//Perform S*z and store the solution in z
	//z now holds St*Minv*S=zaa
	
    Mat66Mult(S,Minv,Minv);	//Perform S*Minv and store the solution in Minv
	//Minv now holds S*Minv, this is done becuase this values is needed later
}

//zab:
//	Function used to calculate z12 and z21 during initialization
//		r is the necessary position vector (r01 for z21 and r02 for z12)
//		z is an empty matrix in which the solution will be stored
//		Minv is the necessary shifter matrix multiplied by the inverse mass matrix
//		S is an empty matrix to be made into the shifter matrix corresponding to the r vector
void zab(float r[],float z[6][6],float Minv[6][6],float S[6][6])
{    
   makeS(S,r);	//Make the shifter matrix

	//The following loop ensures that values that were not set by makeS are set to 0
   for(int c =0; c< 3; c++)
	{
		for(int r=3; r<6; r++)
		{
     			S[r][c]=0;
		}
	}

    Mat66Mult(Minv,S,z);	//Perform Minv*S and store the solution in z
	//z now holds Minv*S=zab
}

//makeS:
//	Function used to create a shifter matrix with the given position vector r
//		S is the matrix in which the shifter matrix will be stored
//		r is the position vector correspoding to the desired shifter matrix
void makeS(float S[6][6],float r[])
{
	//Set the top right 3x3 matrix in S to the skew matrix of r
    S[0][4]=-1*r[2];
    S[0][5]=r[1];
    S[1][3]=r[2];
    S[1][5]=-1*r[0];
    S[2][3]=-1*r[1];
    S[2][4]=r[0];    
}

//MSsetup:
//	Function used to setup the inverse mass matrix and set the shifter matrix to the identity matrix.
//	This function takes advantage of the 2D nature of this simulation by assuming that 
//	the mass matrix will always be diagonal.
//		Minv is the matrix in which the inverse mass matrix will be stored
//		S is the matrix that must be set to the identity matrix to allow for population later
//		m is the mass of the body
//		I is the inertia tensor of the body
void MSsetup(float Minv[6][6],float S[6][6],float m, float I[3][3])
{
   for(int c =0; c< 6; c++)	//Loop through 6 columns
	{
		for(int r=0; r<6; r++)	//Loop through 6 rows
		{
			S[r][c]=0;	//Set the S matrix to 0

			//If the row is less than 3, and the current iteration is on the diagonal of the matrix,
			//then insert the inverse of the inertia tensor as it is needed in M inverse
			if (r<3 && c==r)
			{
    			Minv[r][c]=(1/(I[r][c]));
			}
			//If the row is greater than 2, and the current iteration is on the diagonal of the 
			//of the matrix, then insert 1/m as it is needed in M inverse
			else if(r>2 && c==r)
			{
				Minv[r][c]=(1/m);
			}
			//If not on the diagonal, set that element in M inverse to 0
			else
			{
				Minv[r][c]=0;
			}
			//If on the diagonal, set that element in S to 1
			if(c==r)
			{
				S[r][c]=1;
			}
		}
	}

}

//Fa:
//	Function used to find the state dependant force vector Fa
//		S is a matrix whose first column will hold Fa
//		w is the angular velocity of the body
//		r is the corresponding position vector (r01 or r02)
//		m is the mass of the body
void Fa(float S[6][6],float w, float r[], float m)
{
 	float g = 9.81;	//Set the gravitational constant
	
	//Loop through the first column of S and set it to 0
    for(int rr=0; rr<6; rr++)
	{
         S[rr][0]=0;
    }
	
	//Create Fa based on the force of gravity and the centripedal force on the body
    S[3][0]=-m*w*w*r[0];
    S[4][0]=-g*m-m*w*w*r[1];
}

//r02:
//	Function used to find the r02 vector of the body.  This function takes advantage of the
//	2D nature of this simulation by removing the need for a direction cosine matrix
//		r is the vector in which r02 will be stored
//		angle is the ange of the body with respect to the vertical
//		L is the length of the body
void r02(float r[], float angle,float L)
{
    r[0]=L*sin(angle)/2;
    r[1]=-1*L*cos(angle)/2;
    r[2]=0;

}

//r01:
//	Funtion used to find the r01 position vector of the body. This function takes advantage of the
//	2D nature of this simulation by removing the need for a direction cosine matrix
//		r is the vector in which r02 will be stored
//		angle is the ange of the body with respect to the vertical
void r01(float r[], float angle,float L)
{
    r[0]=-1*L*sin(angle)/2;
    r[1]=L*cos(angle)/2;
    r[2]=0;
}

//transpose:
//	Function used to find the transpose of a 6x6 shifter matrix
//		S is the shifter matrix to be transposed.  The transposed matrix will also be stored here.
void transpose(float S[6][6])
{
	float S2[6][6];	//Temporary matrix used to store values
	
	//Loop through every row and column and save S transpose into S2
	for(int r =0; r< 6; r++)
	{
		for(int c=0; c<6; c++)
		{
			S2[r][c]=S[c][r];
		}
	}

	//Loop through every row and column and save S transpose into S
	for(int r =0; r< 6; r++)
	{
		for(int c=0; c<6; c++)
		{
			S[r][c]=S2[r][c];
		}
	}
}

//Mat66Mult:
//	Function used to multipy two 6x6 matrices.  Note that a temporary 
//	matrix is used to first store the solution. This allows for this 
//	function to store the solution into one of the matrices being multiplied.
//	The form of the operation is A*B=C
//		A is the left matrix
//		B is the right matrix
//		C is the matrix in which the solution is stored
void Mat66Mult(float A[6][6], float B[6][6], float C[6][6])
{
	float j[6][6];	//Temporary matrix used to store the solution

   	for(int r =0; r< 6; r++)	//Loop through all 6 rows
	{
		for(int c=0; c<6; c++)	//Loop through all 6 columns
		{	
			j[r][c]=0;	//Initialize the current element in the temporary matrix to 0

			//The loop below cycles through the row of the first matrix, and the column of the
			//second corresponding to the current row and column of the temporary matrix.  It 
			//takes the dot product of the row and column and stores it in j.
			for(int l =0; l<6; l++)
			{	
				j[r][c]=j[r][c]+(A[r][l]*B[l][c]);
			}
		}
	}

	//Loop through the rows and columns of the temporary matrix, saving the solution into C
	for(int r=0; r<6; r++)
	{
		for(int c=0; c<6; c++)
		{
			C[r][c]=j[r][c];
		}
	}
}

//Mat61Mult:
//	Function used to multipy a 6x6 and a 6x1 matrix.  Note that a temporary 
//	matrix is used to first store the solution. This allows for this function 
//	to store the solution into one of the matrices being multiplied.  It should
//	also be noted that this function takes only 6x6 matrices, and the 6x1 matrix 
//	to be multiplied is assumed to be held in the first column of matrix B.
//	The form of the operation is A*B=C
//		A is the 6x6 matrix
//		B is a 6x6 matrix with the 6x1 matrix stored in its first column
//		C is the matrix in which the solution is stored.  It is stored only in the first column
void Mat61Mult(float A[6][6], float B[6][6], float C[6][6])
{
	float j[6];	//Temporary matrix
 
	for(int r =0; r< 6; r++)	//Loop through all 6 rows
	{	
		j[r]=0;	//Initialize this element of j to 0
		
		//The loop below cycles through the row of the first matrix, and the column of the
		//second corresponding to the current row of the temporary matrix.  It 
		//takes the dot product of the row and column and stores it in j.
		for(int c=0; c<6; c++)
		{
			j[r]+=A[r][c]*B[c][0];
		}
	}	
	
	//Loop through the rows of the temporary matrix, saving the solution into C		
	for(int r = 0; r< 6; r++)
	{
		C[r][0]=j[r];
	}
        
}

