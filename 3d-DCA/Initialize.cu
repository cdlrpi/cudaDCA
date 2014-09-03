//	Initialize.cu
//
//This file contains the function used to find the initial zeta values of every body in the system
#define z11 0
#define z12 6
#define z13 12
#define z21 18
#define z22 24
#define z23 30

#define nz11 0
#define nz12 6
#define nz13 12
#define nz21 13
#define nz22 19
#define nz23 25
//Included Files

#include <stdio.h>
#include <iostream>

//Function Prototypes
//	Functions found in DCAfuncts.cu
void Mult2(double A[6][6], double B[], double C[6][6], int body, int n, int zeta);
void Mult(double A[6][6], double B[6][6], int zeta, int body, int n, double Zetas[],int len);
void S10_Minv(double Mat[6][6], double Mass[], double Inertia[], int n, int body, double Body_Vectors[]);
void get_Minv(double Mat[6][6], double Mass[], double Inertia[], int n, int body, double Body_Vectors[]);
void S20_Minv(double Mat[6][6], double Mass[], double Inertia[], int n, int body, double Body_Vectors[]);
void get_Rots(double Rot1[6][6], double Rot2[6][6], double Angles[], int body);
void Rotate(double Rot1[6][6],double Rot2[6][6], double Zetas[], double Init_Zetas[],int body, int n,double Forces[6][2],bool gravity);
void get_S01( double Body_Vectors[], double S01[6][6], int body);
void get_S02( double Body_Vectors[], double S01[6][6], int body);
void printZeta(double[],int,int,int,int);
void print66(double[6][6]);
void save_DCM(double DCMs[],double temp[3][3],int i,int n);


//Initialize:
//	Function used to find the initial zeta values for every body in the system. n2 is up.
//		state is an array of the state of the system at that timestep
//		oldbds is the list of bodies that do not yet have zeta values
//		newbds is the list of bodies where the new zeta values will be saved
//		n is the number of bodies
void Initialize(double Mass[],double Inertia[], double Zetas[], int n, bool Active_DOF[], double Body_Placement[],double Body_Vectors[],int dof_index[])
{
	double tempS01[6][6];
	double tempS02[6][6];
	double Mat[6][6];
	int count = 0;
	
	for(int body = 0; body<n; body++)
	{	
		dof_index[body]=count;
		for(int idx =0; idx<6; idx++)
		{
			if(Active_DOF[body*6+idx])
			{
				count++;
			}
		}
		get_S01(Body_Vectors,tempS01,body);
		get_S02(Body_Vectors,tempS02,body);
		S10_Minv(Mat,Mass,Inertia,n,body,Body_Vectors);
		print66(tempS01);
		print66(Mat);
		Mult(Mat,tempS01,z11,body,n,Zetas,30);
		printZeta(Zetas,n,body,30,z11);
		Mult(Mat,tempS02,z12,body,n,Zetas,30);
		for(int row=0; row<6; row++)
		{
			for(int col=0; col<6; col++)
			{	
				Zetas[body*30+row*n*30+col+z13]=Mat[row][col];
			}
		}
		S20_Minv(Mat,Mass,Inertia,n,body,Body_Vectors);
		Mult(Mat,tempS01,z21,body,n,Zetas,30);
		Mult(Mat,tempS02,z22,body,n,Zetas,30);

		for(int row=0; row<6; row++)
		{
			for(int col=0; col<6; col++)
			{	
				Zetas[body*30+row*n*30+col+z23]=Mat[row][col];
			}
		}		
	}
}

void update(double Mass[], double Inertia[], double Init_Zetas[], double Zetas[], int n, bool Active_DOF[], double Body_Vectors[], double Coords[], double Speeds[],int dof_index[],double DCM_Angles[])
{

	double Rot1[6][6];
	double Rot2[6][6];
	double Forces[6][2];
		
	for(int body=0; body<n; body++)
	{
		get_Rots(Rot1,Rot2,DCM_Angles,body);
		//get_Force(Forces,Coords,speeds,dof_index,Active_DOF,1);
		Rotate(Rot1,Rot2,Zetas,Init_Zetas,body, n,Forces,1);
	}

}
void serial_operations(bool Active_DOF[], double Coords[], double DCMs[], double speeds[], double omegas[], int dof_index[],int n)
{
	int c;
	double t1 =0;
	double t2 =0;
	double t3 =0;
	double tdot1 =0;
	double tdot2 =0;
	double tdot3 =0;
	double w1 =0;
	double w2 =0;
	double w3 =0;
	for(int i=0; i<n; i++)
	{
		
		c =0;
		if(Active_DOF[i*6])
		{
			t1=Coords[dof_index[i]];
			tdot1=speeds[dof_index[i]];
			c++;
		}
		if(Active_DOF[i*6+1])
		{
			t2=Coords[dof_index[i]+c];
			tdot2=speeds[dof_index[i]+c];
			c++;
		}
		if(Active_DOF[i*6+2])
		{
			t3=Coords[dof_index[i]+c];
			tdot3=speeds[dof_index[i]+c];
		}

		temp[0][0]=cos(t2)*cos(t3);
		temp[0][1]=cos(t1)*sin(t3)+sin(t1)*sin(t2)*cos(t3);
		temp[0][2]=sin(t3)*sin(t1)-cos(t1)*sin(t2)*cos(t3);
		temp[1][0]=-cos(t2)*sin(t3);
		temp[1][1]=cos(t1)*cos(t3)-sin(t1)*sin(t2)*sin(t3);
		temp[1][2]=sin(t1)*cos(t3)+cos(t1)*sin(t2)*sin(t3);
		temp[2][0]=sin(t2);
		temp[2][1]=-sin(t1)*cos(t2);
		temp[2][2]=cos(t1)*cos(t2);
		if(i ==0)
		{
			DCMs[i*3]=temp[0][0];
			DCMs[i*3+1]=temp[0][1];
			DCMs[i*3+2]=temp[0][2];
			DCMs[i*3+n*3]=temp[1][0];
			DCMs[i*3+1+n*3]=temp[1][1];
			DCMs[i*3+2+n*3]=temp[1][2];
			DCMs[i*3+n*6]=temp[2][0];
			DCMs[i*3+1+n*6]=temp[2][1];
			DCMs[i*3+2+n*6]=temp[2][2];
		}
		else
		{
			save_DCM(DCMs,temp,i,n);
		}
		w1+=(tdot1*cos(t2)*cos(t3)+tdot2*sin(t3));
		w2+=(-tdot1*cos(t2)*sin(t3)+tdot2*cos(t3));
		w3+=(tdot3+tdot1*sin(t2);
		omegas[i*3]=w1;
		omegas[i*3+1]=w2;
		omegas[i*3+2]=w3;
}
	

		
			
	}
}

void save_DCM(double DCMs[],double temp[3][3],int i,int n)
{

	
	for(int r =0; r<3; r++)
	{
		for(int c=0; c<3; c++)
		{
			DCMs[i*3+c+r*n*3]=0;
			for(int l =0; l<6; l++)
			{
				DCMs[i*3+c+r*n*3]+=DCMs[(i-1)*3+l+r*n*3]*temp[l][c];
			}
		}
	}
}
			
void get_Force(double Forces[6][2] ,double omegas[] ,double Body_Vectors[], int body);
{
	r1 = Body_Vectors[body*6];
	r2 = Body_Vectors[body*6+1];
	r3 = Body_Vectors[body*6+2];
	w1 = omegas[body*3];
	w2 = omegas[body*3+1];
	w3 = omegas[body*3+2];

	Forces[3][0]= - w2*(r1*w2 - r2*w1) - w3*(r1*w3 - r3*w1);
	Forces[4][0]= w1*(r1*w2 - r2*w1) - w3*(r2*w3 - r3*w2);
	Forces[5][0]= w1*(r1*w3 - r3*w1) + w2*(r2*w3 - r3*w2);
	
	r1 = Body_Vectors[body*6+3];
	r2 = Body_Vectors[body*6+4];
	r3 = Body_Vectors[body*6+5];
	
	Forces[3][1]= - w2*(r1*w2 - r2*w1) - w3*(r1*w3 - r3*w1);
	Forces[4][1]= w1*(r1*w2 - r2*w1) - w3*(r2*w3 - r3*w2);
	Forces[5][1]= w1*(r1*w3 - r3*w1) + w2*(r2*w3 - r3*w2);



}

void get_Rots(double Rot1[6][6], double Rot2[6][6], double Angles[], int body)
{
	double t1=Angles[body*3];
	double t2=Angles[body*3+1];
	double t3=Angles[body*3+2];
	
	Rot1[3][3]=cos(t2)*cos(t3);
	Rot1[3][4]=cos(t1)*sin(t3)+sin(t1)*sin(t2)*cos(t3);
	Rot1[3][5]=sin(t3)*sin(t1)-cos(t1)*sin(t2)*cos(t3);
	Rot1[4][3]=-cos(t2)*sin(t3);
	Rot1[4][4]=cos(t1)*cos(t3)-sin(t1)*sin(t2)*sin(t3);
	Rot1[4][5]=sin(t1)*cos(t3)+cos(t1)*sin(t2)*sin(t3);
	Rot1[5][3]=sin(t2);
	Rot1[5][4]=-sin(t1)*cos(t2);
	Rot1[5][5]=cos(t1)*cos(t2);
	for(int row =0; row<6; row++)
	{
		for(int col=0; col<6; col++)
		{
			Rot2[row][col]=Rot1[col][row];
		}
	}
}
void Rotate(double Rot1[6][6],double Rot2[6][6], double Zetas[], double Init_Zetas[],int body, int n,double Forces[6][2],bool gravity)
{
	double temp[6][6];
	float g = 9.81;

	Mult2(Rot1,	Init_Zetas,temp,body,n, z11);
	Mult(temp,Rot2,nz11	, body, n , Zetas, 26);

	Mult2(Rot1,Init_Zetas,temp,body,n,z12);
	Mult(temp,Rot2,nz12,body,n,Zetas,26);

	Mult2(Rot1,Init_Zetas,temp,body,n,z21);
	Mult(temp,Rot2,nz21,body,n,Zetas,26);

	Mult2(Rot1,Init_Zetas,temp,body,n,z22);
	Mult(temp,Rot2,nz22, body,n,Zetas,26);

	
}
		
void get_S01(double Body_Vectors[],double S01[6][6], int body)
{
	int r1= body*6;
	int r2= r1+1;
	int r3= r2+1;
	for(int row=0; row<6; row++)
	{
		for(int col=0; col<6; col++)
		{
			if(row==col)
			{
				S01[row][col]=1;
			}
			else
			{
			S01[row][col]=0;
			}
		}
	}
	S01[0][4]=-Body_Vectors[r3];
	S01[0][5]=Body_Vectors[r2];
	S01[1][3]=Body_Vectors[r3];
	S01[1][5]=-Body_Vectors[r1];
	S01[2][3]=-Body_Vectors[r2];
	S01[2][4]=Body_Vectors[r1];
}
void get_S02(double Body_Vectors[],double S01[6][6], int body)
{
	int r1= body*6+3;
	int r2= r1+1;
	int r3= r2+1;
	for(int row=0; row<6; row++)
	{
		for(int col=0; col<6; col++)
		{
			if(row==col)
			{
				S01[row][col]=1;
			}
			else
			{
			S01[row][col]=0;
			}
		}
	}
	S01[0][4]=-Body_Vectors[r3];
	S01[0][5]=Body_Vectors[r2];
	S01[1][3]=Body_Vectors[r3];
	S01[1][5]=-Body_Vectors[r1];
	S01[2][3]=-Body_Vectors[r2];
	S01[2][4]=Body_Vectors[r1];
}




void S20_Minv(double Mat[6][6], double Mass[], double Inertia[], int n, int body, double Body_Vectors[])
{
	int inertiaid = body*3;
	int Vectorid = body*6+3;
	Mat[0][0] = 1/Inertia[inertiaid];
	Mat[0][1]=0;
	Mat[0][2]=0;
	Mat[0][3]=0;
	Mat[0][4]=0;
	Mat[0][5]=0;
	Mat[1][0]=0;
	Mat[1][1]=1/Inertia[inertiaid+1];
	Mat[1][2]=0;
	Mat[1][3]=0;
	Mat[1][4]=0;
	Mat[1][5]=0;
	Mat[2][0]=0;
	Mat[2][1]=0;
	Mat[2][2]=1/Inertia[inertiaid+2];
	Mat[2][3]=-0;
	Mat[2][4]=0;
	Mat[2][5]=0;
	Mat[3][0]=0;
	Mat[3][1]=Body_Vectors[Vectorid+2]/Inertia[inertiaid+1];
	Mat[3][2]=-Body_Vectors[Vectorid+1]/Inertia[inertiaid+2];
	Mat[3][3]=(1/Mass[body]);
	Mat[3][4]=0;
	Mat[3][5]=0;
	Mat[4][0]=-Body_Vectors[Vectorid+2]/Inertia[inertiaid];
	Mat[4][1]=0;
	Mat[4][2]=Body_Vectors[Vectorid]/Inertia[inertiaid+2];
	Mat[4][3]=0;
	Mat[4][4]=(1/Mass[body]);
	Mat[4][5]=0;
	Mat[5][0]= Body_Vectors[Vectorid+1]/Inertia[inertiaid];
	Mat[5][1]=-Body_Vectors[Vectorid]/Inertia[inertiaid+1];
	Mat[5][2]=0;
	Mat[5][3]=0;
	Mat[5][4]=0;
	Mat[5][5]=(1/Mass[body]);
}
void get_Minv(double Mat[6][6], double Mass[], double Inertia[], int n, int body, double Body_Vectors[])
{
	int inertiaid = body*3;
	int Vectorid = body*6+3;
	Mat[0][0] = 1/Inertia[inertiaid];
	Mat[0][1]=0;
	Mat[0][2]=0;
	Mat[0][3]=0;
	Mat[0][4]=0;
	Mat[0][5]=0;
	Mat[1][0]=0;
	Mat[1][1]=1/Inertia[inertiaid+1];
	Mat[1][2]=0;
	Mat[1][3]=0;
	Mat[1][4]=0;
	Mat[1][5]=0;
	Mat[2][0]=0;
	Mat[2][1]=0;
	Mat[2][2]=1/Inertia[inertiaid+2];
	Mat[2][3]=-0;
	Mat[2][4]=0;
	Mat[2][5]=0;
	Mat[3][0]=0;
	Mat[3][1]=0;
	Mat[3][2]=0;
	Mat[3][3]=(1/Mass[body]);
	Mat[3][4]=0;
	Mat[3][5]=0;
	Mat[4][0]=0;
	Mat[4][1]=0;
	Mat[4][2]=0;
	Mat[4][3]=0;
	Mat[4][4]=(1/Mass[body]);
	Mat[4][5]=0;
	Mat[5][0]=0;
	Mat[5][1]=0;
	Mat[5][2]=0;
	Mat[5][3]=0;
	Mat[5][4]=0;
	Mat[5][5]=(1/Mass[body]);
}
void S10_Minv(double Mat[6][6], double Mass[], double Inertia[], int n, int body, double Body_Vectors[])
{
	int inertiaid = body*3;
	int Vectorid = body*6;
	Mat[0][0] = 1/Inertia[inertiaid];
	Mat[0][1]=0;
	Mat[0][2]=0;
	Mat[0][3]=0;
	Mat[0][4]=0;
	Mat[0][5]=0;
	Mat[1][0]=0;
	Mat[1][1]=1/Inertia[inertiaid+1];
	Mat[1][2]=0;
	Mat[1][3]=0;
	Mat[1][4]=0;
	Mat[1][5]=0;
	Mat[2][0]=0;
	Mat[2][1]=0;
	Mat[2][2]=1/Inertia[inertiaid+2];
	Mat[2][3]=-0;
	Mat[2][4]=0;
	Mat[2][5]=0;
	Mat[3][0]=0;
	Mat[3][1]=Body_Vectors[Vectorid+2]/Inertia[inertiaid+1];
	Mat[3][2]=-Body_Vectors[Vectorid+1]/Inertia[inertiaid+2];
	Mat[3][3]=(1/Mass[body]);
	Mat[3][4]=0;
	Mat[3][5]=0;
	Mat[4][0]=-Body_Vectors[Vectorid+2]/Inertia[inertiaid];
	Mat[4][1]=0;
	Mat[4][2]=Body_Vectors[Vectorid]/Inertia[inertiaid+2];
	Mat[4][3]=0;
	Mat[4][4]=(1/Mass[body]);
	Mat[4][5]=0;
	Mat[5][0]= Body_Vectors[Vectorid+1]/Inertia[inertiaid];
	Mat[5][1]=-Body_Vectors[Vectorid]/Inertia[inertiaid+1];
	Mat[5][2]=0;
	Mat[5][3]=0;
	Mat[5][4]=0;
	Mat[5][5]=(1/Mass[body]);
}


void Mult(double A[6][6], double B[6][6], int zeta, int body, int n, double Zetas[],int len)
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
			Zetas[zeta+body*len+r*n*len+c]=j[r][c];
		}
	}
}
			
void Mult2(double A[6][6], double B[], double C[6][6], int body, int n, int zeta)
{
   	for(int r =0; r< 6; r++)	//Loop through all 6 rows
	{
		for(int c=0; c<6; c++)	//Loop through all 6 columns
		{	
			C[r][c]=0;	//Initialize the current element in the temporary matrix to 0

			//The loop below cycles through the row of the first matrix, and the column of the
			//second corresponding to the current row and column of the temporary matrix.  It 
			//takes the dot product of the row and column and stores it in j.
			for(int l =0; l<6; l++)
			{	
				C[r][c]=C[r][c]+(A[r][l]*B[l*30*n+30*body+c+zeta]);
			}
		}
	}

}								
void print66(double A[6][6])
{
	std::cout<<"\n\n\n";
	for(int r =0; r<6; r++)
	{
		for(int c =0; c<6; c++)
		{
			std::cout<<A[r][c]<<'\t';
		}
		std::cout<<std::endl;
	}
}	
