#include "classes.h"
#include <iostream>
void CudaInitialize(InitBody* oldbs, float x[], int n, Body *newbs);
void printm(float A[6][6]);
void printa(float A[], int x);
void cudaAssemble(Body *bodies, int num, Body *newbds, int odd, int newlen);
void RecDCA(Body *bodies, int n, int i, float AF[]);
void cudaDisassemble(float OldAF[], Body *morebds, Body *lessbds, int odd, int morelen, int lesslen, float AF[]);
void solve_BCs(Body *bodies, float AF[]);
void printA(float M[6][10]);
//Function that helps DCA
void printA(float M[6][10])
{
	for(int i = 0; i < 6; i++)
	{
		for( int j = 0; j<10;j++)
		{
			std::cout<<M[i][j]<<"\t";
		}
		std::cout<<std::endl;
	}
}
void DCAhelp(float state[], InitBody *bs, Joint *js,int n, float Y[])
{
	float *AF = (float*) malloc(sizeof(float)*n*4*6);
	float A[6][n*2];
	Body *bodies = new Body[n];
	CudaInitialize(bs, state, n, bodies);


	RecDCA(bodies, n, 0, AF);

	for(int c = 0, c2=0 ; c<2*n; c++, c2+=2)
	{
		for(int r = 0; r<6;r++)
		{
			A[r][c]=AF[c2+r*4*n];
		}
	}
	Y[n]=A[2][0];

	for(int i = n+1, j=2; i<n*2; i++, j+=2)
	{
		Y[i]= A[2][j]-A[2][j-1];	
	}
	for(int i = 0; i<n; i++)
	{
		Y[i]=state[i+n];
	}
	free(AF);
	delete[] bodies;	
}

void RecDCA(Body *bodies, int n, int i, float AF[])
{
	if (n==1)
	{
		solve_BCs(bodies, AF);
		
		
	}
	else
	{	
		
		int newlen;
		
		int odd = 0;
	
		Body *newbds;
		if(n % 2 == 0)
		{
			newlen = (int) (n/2);
		
		}
		else
		{
			newlen = (int)((n+1)/2);
			odd = 1;
		}
		
		newbds = new Body[newlen];
		cudaAssemble(bodies, n, newbds, odd, newlen);
		float *AFo=(float*)malloc(sizeof(float)*6*newlen*4);
		RecDCA(newbds,newlen,i+1 ,AFo);
		cudaDisassemble(AFo, bodies , newbds, odd, n, newlen, AF);
		delete[] newbds;
		free(AFo); 
		
	}
}
		
			

