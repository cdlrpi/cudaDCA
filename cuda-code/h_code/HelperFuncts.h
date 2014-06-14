#include "classes.h"
void printm(float A[6][6])
{	
	for(int r = 0;r<6;r++)
	{
		for(int c = 0; c < 6; c++)
		{

			std::cout<<A[r][c]<<"\t";
		}
		std::cout<<std::endl;
	}
	std::cout<<std::endl<<std::endl;
}
void printa(float A[], int x)
{
	for(int i = 0; i<x;i++)
	{
		std::cout<<A[i]<<std::endl;
	}
	std::cout<<std::endl<<std::endl;
}
void matmult61(float A[6][6], float B[6], float C[6])
{	
	float j;
	for(int r = 0; r<6; r++)
	{
		j=0;
		for(int c = 0; r<6; r++)
		{
			j+= (A[r][c]*B[c]);
		}
		C[r] = j;
	}
}
	
void pend_init(InitBody *bs,int n,float mass, float length)
{
	int k=0;
	float II;
	while(k<n)
	{
		bs[k].m=mass;
		bs[k].l=length;
		II= mass*(pow(length,2)/12.0);
		bs[k].I[0][0]=II;
		bs[k].I[1][1]=1;
		bs[k].I[2][2]=II;
		k++;
	}
}

//Function to set the initial conditions to that of a full drop from 90 degrees.
void full_drop(float x[],int n)
{
	int k = 0;
	x[k]=3.14159/2;
	k++;
	while(k<2*n)
	{
		x[k]=0;
		k++;
	}
}

//Function to initialize the time matrix
void Time_init(float time[], float step, float final)
{
	float t = 0;
	int c = 0;
	while (t<final)
	{
		time[c]=t;
		t+=step;
		c++;
	}
}

//Function to copy array B into array A
void arycpy(float A[],float B[],int n)
{
	int i=0;
	while (i<n)
	{
		A[i]=B[i];
		i++;
	}
}

//Function to copy array B into the second half of array A
void arycpy2(float A[], float B[], int n)
{
	int i = 0;
	while(i<n)
	{
		A[i+n]=B[i];
		i++;
	}
}

//Function to copy the second half of array B into array A
void arycpy3(float A[], float B[], int n)
{
	int i = 0;
	while(i<n)
	{
		A[i]=B[i+n];
		i++;
	}
}

//Function to set up the qs and qdots
void set_up(float A[], float B[], float C[], int n , float h)
{
	int i = 0;
	while (i<n)
	{
		C[i]=A[i]+(h*B[i]);
	}
}

//Function to get the final q to be used in the next timestep
void get_final_q(float x[], float h, float a[], float b[], float c[], float d[], float R[], int n)
{
	int i = 0;
	while(i<n)
	{
		R[i]= x[i] + ((h/3.0)*(a[i]+(2*b[i])+(2*c[i])+d[i]));
	}
}

//Function to get the final qdot to be used in the next timestep
void get_final_qdot(float x[], float h, float a[], float b[], float c[], float d[], float R[], int n)
{
	int i = 0;
	while(i<n)
	{
		R[i+n]= x[i] + ((h/3.0)*(a[i]+(2*b[i])+(2*c[i])+d[i]));
	}
}

//Joint initialization function
//The joint is assumed to be that of a 2d pendulum
Joint::Joint()
{
	//int P[]={0,0,1,0,0,0};
	//int D[6][5]={{1,0,0,0,0},{0,1,0,0,0},{0,0,0,0,0},{0,0,1,0,0},{0,0,0,1,0},{0,0,0,0,1}};
}
