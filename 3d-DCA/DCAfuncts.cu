//	DCAfuncts.cu
//
//This file contains all of the functions that are used during initialization, 
//assembly, and disassembly.

#include <iostream>
//Function Prototypes
//	All function called in this file are written in this file
void transpose(double S[6][6]);
void Mat66Mult(double A[6][6], double B[6][6], double C[6][6]);
void Mat61Mult(double A[6][6], double B[6][6], double C[6][6]);
void makeS(double S[6][6],double r[]);
void set_D(double D[6][6]);
void set_Dt(double D[6][6]);
void get_D2(int Pindex[], int numbods, int b, double D[],double A[6][6]);
void printm(double A[6][6]);
void rotate_PD(double DCMs[], double D[6][6], int n, int b);
void rotate_D(double DCMs[], double Ds[], int n, int b);
void makePdotU(double PdotUs[], int n, int b, double P[6][6],bool active_DOF[], double ws[], double speeds[], int dof_index[]);
//get_X:
//	Function used to find the intermediate quantity X for a body
//	Considering the two bodies used to create X as b1 and b2:
//		z1 is a matrix containing z11 for b2
//		z2 is z22 for b1
//		D is the matrix in which X will be stored
void get_X(double z1[6][6], double z2 [6][6], double Ds[],int n, bool active_DOF, int b,int Pindex[],int numbods,double A[6][6])
{	
	for(int r = 0; r<6; r++)	//Loop through all 6 rows
	{
		for(int c=0; c<6; c++)	//Loop through all 6 columns
		{
    		z1[r][c]=z1[r][c]+z2[r][c];	//Add z1 and z2 together and store the solution in z1
		}
	}
	//Completing the above loop, z1 now holds b1.z22+b2.z11
    get_D2(Pindex, numbods, b,Ds,A);	//put the D matrix corresponding to a joint in this simulation in D

    Mat66Mult(z1,A,z1);	//Perform z1*D and store the solution in z1
	//z1 now holds (b1.z22+b2.z11)*D

    transpose(A);	//put the transpose of the D matrix in D
    Mat66Mult(A,z1,A); //Perform D*z1 and store the solution in D

	//D now holds Dt*(b1.z22+b2.z11)*D=X
}
void get_D2(int Pindex[], int numbods, int b, double Ds[],double D[6][6])
{
	for(int r =0; r<6; r++)
	{
		for(int c =0; c<6; c++)
		{
			D[r][c]=Ds[Pindex[b]*6+c+r*6*numbods];
		}
	}
}
void get_P(double Ps[], double PdotUs[], int n, bool active_DOF[], double DCMs[], int b,double ws[], double speeds[], int dof_index[])
{
	double P[6][6];
	int c =0;
	int r =0;
	for(int x =0; x<6; x++)
	{
		for(int y=0; y<6; y++)
		{
			P[x][y]=0;
		}
	}
	while(r<6)	//Loop through all 6 rows
	{
    	if(active_DOF[b*6+r])
		{
			P[r][c]=1;
			c++;			
		}
		r++;				
	}

	rotate_PD(DCMs, P, n, b);
	for(int r =0; r<6; r++)
	{
		for(int c =0; c<6; c++)
		{
			Ps[c+b*6+n*r*6]=P[r][c];
		}
	}
	makePdotU(PdotUs,n,b,P,active_DOF,ws,speeds,dof_index);
}


void makePdotU(double PdotUs[], int n, int b, double P[6][6],bool active_DOF[], double ws[], double speeds[], int dof_index[])
{
	double U[6];
	double temp[6][6];
	int x=0;
	for(int i =0; i<6; i++)
	{
		U[i]=0;
		if(active_DOF[b*6+i])
		{
			U[x]=speeds[dof_index[b]+x];
			x++;
		}
		for(int j =0; j<6; j++)
		{
			temp[i][j]=0;
		}
	}
	for(int i =0; i<6; i++)
	{
		for(int j=0; j<4; j+=3)
		{
			temp[0+j][i]=(ws[b*3+1]*P[2+j][i]-ws[b*3+2]*P[1+j][i]);
			temp[1+j][i]=-1*(ws[b*3]*P[2+j][i]-ws[b*3+2]*P[0+j][i]);
			temp[2+j][i]=(ws[b*3]*P[1+j][i]-ws[b*3+1]*P[0+j][i]);
		}
	}
		//Temporary matrix
 
	for(int r =0; r< 6; r++)	//Loop through all 6 rows
	{	
		PdotUs[b*6+r]=0;	//Initialize this element of j to 0
		
		//The loop below cycles through the row of the first matrix, and the column of the
		//second corresponding to the current row of the temporary matrix.  It 
		//takes the dot product of the row and column and stores it in j.
		for(int c=0; c<6; c++)
		{
			PdotUs[b*6+r]+=temp[r][c]*U[c];
		}
		//std::cout<<PdotUs[b*6+r];
	}
}	
		
		

//set_D:
//	Function used to set the matrix that defines the subspace of prohibited motion.
//		D is the matrix in which the prohibited motion matrix will be stored.  
void get_D(double Ds[],int n, bool active_DOF[], double DCMs[], int b)
{
	int c =0;
	int r =0;
	for(int x =0; x<6; x++)
	{
		for(int y=0; y<6; y++)
		{
			Ds[b*6+y+n*6*x]=0;
		}
	}
	while(r<6)	//Loop through all 6 rows
	{
    	if(!active_DOF[b*6+r])
		{
			Ds[r*n*6+c+b*6]=1;
			c++;			
		}
		r++;				
	} 

	rotate_D(DCMs, Ds, n, b);

}
void rotate_D(double DCMs[], double Ds[], int n, int b)
{
	for(int c =0; c<6; c++)
	{
		for(int r=0; r<3; r++)
		{
			if(Ds[c+b*6+r*n*6])
			{
				for(int r2=0; r2<3; r2++)
				{
					Ds[r2*n*6+c+b*6]=DCMs[r+b*3+r2*n*3];
				}
				r=3;
			}	
		}
	}
	for(int c =0; c<6; c++)
	{
		for(int r = 3; r<6; r++)
		{
			if(Ds[r*n*6+c+b*6])
			{
				for(int r2=0; r2<3; r2++)
				{
					Ds[(r2+3)*n*6+c+b*6]=DCMs[r2*n*3+b*3+r-3];
				}
				r=6;
			}
		}
	}
}
void rotate_PD(double DCMs[], double D[6][6], int n, int b)
{
	for(int c =0; c<6; c++)
	{
		for(int r=0; r<3; r++)
		{
			if(D[r][c])
			{
				for(int r2=0; r2<3; r2++)
				{
					D[r2][c]=DCMs[r+b*3+r2*n*3];
				}
				r=3;
			}	
		}
	}
	for(int c =0; c<6; c++)
	{
		for(int r = 3; r<6; r++)
		{
			if(D[r][c])
			{
				for(int r2=0; r2<3; r2++)
				{
					D[r2+3][c]=DCMs[r-3+b*3+r2*n*3];
				}
				r=6;
			}
		}
	}
}
				
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
//set_Dt:
//	Function used to set the transpose of the matrix that defines the subspace of prohibited motion.
//		D is the matrix in which the transposed prohibited motion matrix will be stored.  
/*
void transpose(double D[6][6])
{
	double temp[6][6];
    for(int r = 0; r<6; r++)	//Loop through all 6 rows
	{
		for(int c=0; c<6; c++)	//Loop through all 6 columns
		{
    		temp[c][r] = D[r][c];	//Set the D matrix to 0;
		}
	}
	for(int r =0; r<6; r++)
	{
		for(int c =0; c<6; c++)
		{
			D[r][c]=temp[r][c];
		}
	}	
}
*/
//invert_X:
//	Function used to invert the X matrix.  This function takes advantage of the nature of this 
//	2D simulation by only inverting the 2x2 matrix inside of X that will effect the solution.
//	This was done for efficiency of the simulation as well as ease of implementation in this case.
//		X is the intermediate quantity X; X inverse will be stored in X, replacing it
void invert_X(double X[6][6])
{
	double temp1;
	double temp2;
	double eye[5][5];
	for(int r =0; r<5; r++)
	{
		for(int c =0; c<5; c++)
		{
			eye[r][c]=0;
			if(c==r)
			{
				eye[r][c]=1;
			}
		}
	}
	//Variable Declarations
	for(int r =0; r<5; r++)
	{
		temp1=X[r][r];
		for(int c =0; c<5; c++)
		{
			X[r][c]=X[r][c]/temp1;
			eye[r][c]=eye[r][c]/temp1;
		}
		for(int i=0; i<5; i++)
		{
			if(i!=r)
			{
				temp2=X[i][r];
				for(int c =0; c<5; c++)
				{
					X[i][c]-=(X[r][c]*temp2);
					eye[i][c]-=(eye[r][c]*temp2);
				}
			}
		}
	}
	for(int r =0; r<5; r++)
	{
		for(int c =0; c<5; c++)
		{
			X[r][c]=eye[r][c];
		}
	}
}

//make_W:
//	Function used to create the intermediate quantity W
//		Xinv is Xinverse;  W will be stored here
//		D is an empty matrix used for matrix calculations
void make_W(double Xinv[6][6], double Ds[], int n, int b, int numbods, int Pindex[])                    
{ 
	double D[6][6];  
	get_D2(Pindex, numbods, b,Ds,D);

    //Perform Xinv*D and store the solution in Xinv
	//Xinv now holds Xinv*Dt
	Mat66Mult(D,Xinv,Xinv);
    transpose(D);	//Put the D matrix in D
    	//Perform D*Xinv and store the solution in Xinv
	 Mat66Mult(Xinv,D,Xinv);//Xinv now holds D*Xinv*Dt=W
}

//zaa:
//	Function used to calculate z11 and z22 during initialization
//		r is the necessary position vector (r01 for z11 and r02 for z22)
//		z is the matrix in which the solution will be stored
//		Minv is the inverse mass matrix
//		S is an empty matrix that is made into the shifter matrix corresponding to the r vector
void zaa(double r[],double z[6][6],double Minv[6][6],double S[6][6])
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
void zab(double r[],double z[6][6],double Minv[6][6],double S[6][6])
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
void makeS(double S[6][6],double r[])
{
	//Set the top right 3x3 matrix in S to the skew matrix of r
    S[0][4]=-1*r[2];
    S[0][5]=r[1];
    S[1][3]=r[2];
    S[1][5]=-1*r[0];
    S[2][3]=-1*r[1];
    S[2][4]=r[0];    
}

void makeSt(double S[6][6],double r[])
{
	//Set the top right 3x3 matrix in S to the skew matrix of r
    S[4][0]=-1*r[2];
    S[5][0]=r[1];
    S[3][1]=r[2];
    S[5][1]=-1*r[0];
    S[3][2]=-1*r[1];
    S[4][2]=r[0];    
}

//MSsetup:
//	Function used to setup the inverse mass matrix and set the shifter matrix to the identity matrix.
//	This function takes advantage of the 2D nature of this simulation by assuming that 
//	the mass matrix will always be diagonal.
//		Minv is the matrix in which the inverse mass matrix will be stored
//		S is the matrix that must be set to the identity matrix to allow for population later
//		m is the mass of the body
//		I is the inertia tensor of the body
void MSsetup(double Minv[6][6],double S[6][6],double m, double I[3][3])
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
void Fa(double S[6][6],double w, double r[], double m)
{
 	double g = 9.81;	//Set the gravitational constant
	
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
void r02(double r[],double L)
{
    r[0]=0;
    r[1]=-L/2;
    r[2]=0;

}

//r01:
//	Funtion used to find the r01 position vector of the body. This function takes advantage of the
//	2D nature of this simulation by removing the need for a direction cosine matrix
//		r is the vector in which r02 will be stored
//		angle is the ange of the body with respect to the vertical
void r01(double r[], double L)
{
    r[0]=0;
    r[1]=L/2;
    r[2]=0;
}

//transpose:
//	Function used to find the transpose of a 6x6 shifter matrix
//		S is the shifter matrix to be transposed.  The transposed matrix will also be stored here.
void transpose(double S[6][6])
{
	double S2[6][6];	//Temporary matrix used to store values
	
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
void Mat66Mult(double A[6][6], double B[6][6], double C[6][6])
{
	double j[6][6];	//Temporary matrix used to store the solution

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
void Mat61Mult(double A[6][6], double B[6][6], double C[6][6])
{
	double j[6];	//Temporary matrix
 
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

void Mat61Mult2(double A[6][6], double B[6], double C[6])
{
	double j[6];
	
	for(int r =0; r<6; r++)
	{
		j[r]=0;
		for(int c =0; c<6; c++)
		{
			j[r]+=A[r][c]*B[c];
		}
	}
	for(int r =0; r<6; r++)
	{
		C[r]=j[r];
	}
}

void update(double R[6][6], double Rt[6][6], double zeta1[], double zeta2[], int n, int bodynum, int zetanum)
{
	double temp[6][6];
	for(int r =0; r<6; r++)
	{
		for(int c =0; c<6; c++)
		{
			temp[r][c]=0;
			for(int i =0; i<6; i++)
			{
				temp[r][c]+=zeta1[i+r*26*n+26*bodynum+zetanum]*Rt[i][c];
			}
			//std::cout<<Rt[r][c]<< "  ";
		}
		//std::cout<<std::endl;
	}

	for(int r =0; r<6; r++)
	{
		for(int c =0; c<6; c++)
		{
			zeta2[c+r*26*n+26*bodynum+zetanum]=0;
			for(int i =0; i<6; i++)
			{
				zeta2[c+r*26*n+26*bodynum+zetanum]+=R[r][i]*temp[i][c];
			}
		}
	}
}

void savez3(double Minv[6][6],double S[6][6] ,double Fa[], double nZetas[], int n, int bodynum, int zetanum)
{
	double temp[6][6];
	Mat66Mult(S,Minv,temp);
	for(int r =0; r<6; r++)
	{
		nZetas[r*26*n+bodynum*26+zetanum]=0;
		for(int c =0; c<6; c++)
		{
			nZetas[r*26*n+bodynum*26+zetanum]+=temp[r][c]*Fa[c];
		}
	}
}

void makeR(double R[6][6], double Rt[6][6], double angle)
{
	//std::cout<<angle<<std::endl;
	for(int r =0; r<6; r++)
	{
		for(int c =0; c<6; c++)
		{
			R[r][c]=0;
			Rt[r][c]=0;
		}
		R[r][r]=1;
		Rt[r][r]=1;
	}
	//std::cout<<std::endl<<cos(angle)<<std::endl;
	R[3][3]=cos(angle);
	R[4][4]=cos(angle);
	Rt[3][3]=cos(angle);
	Rt[4][4]=cos(angle);
	R[3][4]=-sin(angle);
	Rt[3][4]=sin(angle);
	R[4][3]=sin(angle);
	Rt[4][3]=-sin(angle);

}


			

