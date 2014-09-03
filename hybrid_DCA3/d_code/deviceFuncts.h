//	deviceFuncts.h
//
//This file contains all of the functions that are used during initialization, 
//assembly, and disassembly.  All functions in this file run on the gpu.

//Function Prototypes
//	All function called in this file are written in this file
__device__ void set_Dt_gpu(double D[6][6],int row,int col);
__device__ void set_D_gpu(double D[6][6],int row,int col);
__device__ void transpose_gpu(double S[6][6],int row, int col);
__device__ void makeS_gpu(double S[6][6],double *r);
__device__ void Mat61Mult_gpu(double A[6][6], double B[6][6], double C[6][6],int row, int col);
__device__ void Mat66Mult_gpu(double A[6][6], double B[6][6], double C[6][6],int row, int col);

//get_X:
//	Function used to find the intermediate quantity X for a body
//	Considering the two bodies used to create X as b1 and b2:
//		z1 is a matrix containing z11 for b2
//		z2 is z22 for b1
//		D is the matrix in which X will be stored
//		row and col are the threadIdx.y and threadIdx.x
__device__ void get_X_gpu(double z1[6][6], double z2 [6][6], double D[6][6],int row , int col)
{
	//Add z1 and z2 together and store the solution in z1
    z1[row][col]=z1[row][col]+z2[row][col];
    
    set_D_gpu(D,row,col);	//put the D matrix corresponding to a joint in this simulation in D

    Mat66Mult_gpu(z1,D,z1,row,col);	//Perform z1*D and store the solution in z1
	//z1 now holds (b1.z22+b2.z11)*D

    set_Dt_gpu(D,row,col);	//put the transpose of the D matrix in D

    Mat66Mult_gpu(D,z1,D,row,col);//Perform D*z1 and store the solution in D
	//D now holds Dt*(b1.z22+b2.z11)*D=X	
}

//set_D:
//	Function used to set the matrix that defines the subspace of prohibited motion.
//		D is the matrix in which the prohibited motion matrix will be stored. 
//		row and col are	the threadIdx.y and threadIdx.x
__device__ void set_D_gpu(double D[6][6],int row,int col)
{
    D[row][col]=0;	//Set the D matrix to 0

    if(row+col==0)	//Ensure only one thread enters the if statement
    {
		//Set the D matrix according to the 2D joints used in this simulation
        D[0][0]=1;
        D[1][1]=1;
        D[3][2]=1;
        D[4][3]=1;
        D[5][4]=1;
    }
	__syncthreads();	//Ensure D is done being made before continuing
}

//set_Dt:
//	Function used to set the transpose of the matrix that defines the subspace of prohibited motion.
//		D is the matrix in which the transposed prohibited motion matrix will be stored.
//		row and col are	the threadIdx.y and threadIdx.x
__device__ void set_Dt_gpu(double D[6][6],int row,int col)
{
    D[row][col]=0;	//Set the D matrix to 0;

    if(row+col==0)	//Ensure only one thread enters the if statement
    {
		//Set the D matrix according to the 2D joints used in this simulation
        D[0][0]=1;
        D[1][1]=1;
        D[2][3]=1;
        D[3][4]=1;
        D[4][5]=1;
    }
	__syncthreads();	//Ensure D is done being made before continuing
}

//invert_X:
//	Function used to invert the X matrix.  
//		X is the intermediate quantity X; X inverse will be stored in X, replacing it
//		A is a temporary matrix used for matrix operations
//		row and col are	the threadIdx.y and threadIdx.x
__device__ void invert_X_gpu(double X[6][6], double A[6][6],int row, int col)
{/*
     int s;
     int j;
     double temp;
     double temp2;
     A[row][col]=0;
     if (row==col)
     {
          A[row][col]=1;
         for(s=0;s<5;s++)
         {
            temp=X[s][s];
            X[s][col]=(X[s][col]/temp);
            A[s][col]=(A[s][col]/temp);
            for(j=0;j<5;j++)
            {
                if(s != j)
                {
                     temp2=X[j][s];
                     X[j][col]=X[j][col]-(X[s][col]*temp2);
                     A[j][col]=A[j][col]-(A[s][col]*temp2);
                }
            }
          }
        }
        __syncthreads();
        X[row][col]=A[row][col];
        __syncthreads();*/
//Variable Declarations
	double det;	//The determinant of the 2x2 matrix
	__shared__ double Xinv[2][2];	//The inverse of the 2x2 matrix

	//Matrix inversion
	if(threadIdx.x==0 && threadIdx.y==0)
	{
	det=(X[2][2]*X[3][3])-(X[2][3]*X[3][2]);	//Find the determinant of the 2x2 matrix

	//Create the inverse matrix based on the trick for 2x2 matrices 
	Xinv[0][0]=X[3][3]/det;
	Xinv[0][1]=-1*X[2][3]/det;
	Xinv[1][0]=-1*X[3][2]/det;
	Xinv[1][1]=X[2][2]/det;
	}
	__syncthreads();
	X[row][col]=0;

//Save the Solution
		if(row<4 && col<4 && row>=2 && col >= 2)
		{
			X[row][col]=Xinv[row-2][col-2];	//Save X inverse in X
		}
	__syncthreads();


}

//make_W:
//	Function used to create the intermediate quantity W
//		Xinv is Xinverse;  W will be stored here
//		D is an empty matrix used for matrix calculations
//		row and col are	the threadIdx.y and threadIdx.x
__device__ void make_W_gpu(double Xinv[6][6], double D[6][6],int row , int col)                    
{
    if (row == 5 || col==5)
    {
         Xinv[row][col]=0;	//Set the values in row and column 6 to 0
    }

    set_Dt_gpu(D,row,col);	//Put the D transpose matrix in D
    Mat66Mult_gpu(Xinv,D,Xinv,row,col);	//Perform Xinv*D and store the solution in Xinv
	//Xinv now holds Xinv*Dt

    set_D_gpu(D,row,col);	//Put the D matrix in D
    Mat66Mult_gpu(D,Xinv,Xinv,row,col);	//Perform D*Xinv and store the solution in Xinv
	//Xinv now holds D*Xinv*Dt=W
}

//zaa:
//	Function used to calculate z11 and z22 during initialization
//		r is the necessary position vector (r01 for z11 and r02 for z22)
//		z is the matrix in which the solution will be stored
//		Minv is the inverse mass matrix
//		S is an empty matrix that is made into the shifter matrix corresponding to the r vector
//		row and col are	the threadIdx.y and threadIdx.x
__device__ void zaa_gpu(double r[],double z[6][6],double Minv[6][6],double S[6][6], int row, int col)
{
     
     if ((row+col)==0)	//Ensure only 1 thread enters the if
     {
          makeS_gpu(S,r);	//Create the shifter matrix
     }

   	__syncthreads();	//Ensure S is done being created
   
     Mat66Mult_gpu(Minv,S,z,row,col);	//Perform Minv*S and store the solution in z
	//z now holds Minv*S
     
     transpose_gpu(S,row,col);	//Transpose the shifter matrix
     Mat66Mult_gpu(S,z,z,row,col);	//Perform S*z and store the solution in z
	//z now holds St*Minv*S=zaa

     Mat66Mult_gpu(S,Minv,Minv,row,col);	//Perform S*Minv and store the solution in Minv
	//Minv now holds S*Minv, this is done becuase this values is needed later
}

//zab:
//	Function used to calculate z12 and z21 during initialization
//		r is the necessary position vector (r01 for z21 and r02 for z12)
//		z is an empty matrix in which the solution will be stored
//		Minv is the necessary shifter matrix multiplied by the inverse mass matrix
//		S is an empty matrix to be made into the shifter matrix corresponding to the r vector
//		row and col are	the threadIdx.y and threadIdx.x
__device__ void zab_gpu(double r[],double z[6][6],double Minv[6][6],double S[6][6],int row , int col)
{    
    if (row+col==0)	//Ensure only 1 thread enters the if
	{
          makeS_gpu(S,r);	//Create the shifter matrix
	}

	if (row>2 && col<3)
	{
		S[row][col]=0;	//Ensure the proper elements in S are 0
	}
	__syncthreads();	//Ensure S is done being made before continuing
    Mat66Mult_gpu(Minv,S,z,row,col);//Perform Minv*S and store the solution in z
	//z now holds Minv*S=zab	
}

//makeS:
//	Function used to create a shifter matrix with the given position vector r
//		S is the matrix in which the shifter matrix will be stored
//		r is the position vector correspoding to the desired shifter matrix
__device__ void makeS_gpu(double S[6][6],double r[])
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
//		row and col are	the threadIdx.y and threadIdx.x
//		Iindex is the index used to read the inertia tensor
__device__ void MSsetup_gpu(double Minv[6][6],double S[6][6],int row , int col,double m, double I[],int Iindex)
{
    Minv[row][col]=0;	//Set Minv to 0
    S[row][col]=0;	//Set S to 0

	//Fill the inverse mass matrix appropriately
    if (row<3 && col==row)
    {
        Minv[row][col]=1/(I[Iindex]);	
        Minv[row+3][col+3]=1/m;
    }

    if (col==row)	//Ensure only 1 thread enters the if
    {
        S[row][col]=1;	//Make S into the identity matrix
    }

	__syncthreads();	//Ensure S and M are done being made before coninuing  
}

//Fa:
//	Function used to find the state dependant force vector Fa
//		S is a matrix whose first column will hold Fa
//		w is the angular velocity of the body
//		r is the corresponding position vector (r01 or r02)
//		m is the mass of the body
//		row and col are	the threadIdx.y and threadIdx.x
__device__ void Fa_gpu(double S[6][6],double w, double r[], double m,int row, int col)
{
	double g = 9.81;	//Set the gravitational constant

	if (col == 0)
	{
		S[row][col]=0;	//set the first column of S to 0
	}

	if (row +col == 0)	//Ensure only 1 thread enters the if
	{
		//Create Fa based on the force of gravity and the centripedal force on the body
		S[3][0]=-m*w*w*r[0];
		S[4][0]=-g*m-m*w*w*r[1];
	}
	__syncthreads();	//Ensure Fa is done being made before continuing
}

//r02:
//	Function used to find the r02 vector of the body.  This function takes advantage of the
//	2D nature of this simulation by removing the need for a direction cosine matrix
//		r is the vector in which r02 will be stored
//		angle is the ange of the body with respect to the vertical
//		L is the length of the body
__device__ void r02_gpu(double r[], double angle,double L)
{
    r[0]=L*sin(angle)/2;
    r[1]=-1*L*cos(angle)/2;
    r[2]=0;

	__syncthreads();
}

//r01:
//	Funtion used to find the r01 position vector of the body. This function takes advantage of the
//	2D nature of this simulation by removing the need for a direction cosine matrix
//		r is the vector in which r02 will be stored
//		angle is the ange of the body with respect to the vertical
__device__ void r01_gpu(double r[], double angle,double L)
{
    r[0]=-1*L*sin(angle)/2;
    r[1]=L*cos(angle)/2;
    r[2]=0;

	__syncthreads();
}

//transpose:
//	Function used to find the transpose of a 6x6 shifter matrix
//		S is the shifter matrix to be transposed.  The transposed matrix will also be stored here.
//		row and col are	the threadIdx.y and threadIdx.x
__device__ void transpose_gpu(double S[6][6], int row, int col)
{
	if (row<3 && col>2)
	{
		S[col][row]=S[row][col];
		S[row][col]=0;
	}

	__syncthreads();
}

//Mat66Mult:
//	Function used to multipy two 6x6 matrices.  Note that a temporary 
//	matrix is used to first store the solution. This allows for this 
//	function to store the solution into one of the matrices being multiplied.
//	The form of the operation is A*B=C
//		A is the left matrix
//		B is the right matrix
//		C is the matrix in which the solution is stored
//		row and col are	the threadIdx.y and threadIdx.x
__device__ void Mat66Mult_gpu(double A[6][6], double B[6][6], double C[6][6],int row, int col)
{
    double j = 0;	//Temporary value used to store the solution
    int l=0;	//counter
    while (l<6)
    {
        j+=A[row][l]*B[l][col];	//dot product of row and column
        l++;
    }
    __syncthreads();
    C[row][col]=j;
    __syncthreads();
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
//		row and col are	the threadIdx.y and threadIdx.x
__device__ void Mat61Mult_gpu(double A[6][6], double B[6][6], double C[6][6],int row, int col)
{
    double j = 0;//Temporary value used to store the solution
    int l=0;	//counter
    if (col==0)
    {
        while (l<6)
        {
            j+=A[row][l]*B[l][0];	//dot product of row and column
            l++;
        }
        C[row][0]=j;
    }
    __syncthreads();
}
/*
//////////////////////////////////////////////////////////////////////////////////////////////////////
//Functions below this line are used for debugging purposes only and can be deleted at any time
__device__ void printit_gpu(double A[6][6],int n)
{
	if(threadIdx.x==threadIdx.y && threadIdx.y==0)
	{
		for(int i =0; i<n; i++)
		{
		if(blockIdx.x==i)
		{
		for(int r = 0; r<6; r++)
		{
			for(int c=0; c<6; c++)
			{
				printf("%f  ",A[r][c]);
			}
			printf("\n");
		}
		printf("\n");
	}
	
}
}}
*/
