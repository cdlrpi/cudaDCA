#include "classes.h"

void CudaInitialize(InitBody* oldbs, float x[], int n, Body *newbs);
void printm(float A[6][6]);
void printa(float A[], int x);
void cudaAssemble(Body *bodies, int num, Body *newbds, int odd, int newlen);
void RecDCA(Body *bodies, int n, int i, float AF[]);
void cudaDisassemble(float OldAF[], Body *morebds, Body *lessbds, int odd, int morelen, int lesslen, float AF[]);
void solve_BCs(Body *bodies, float AF[]);
//Function that helps DCA
void DCAhelp(float state[], InitBody *bs, Joint *js,int n, float Y[])
{
	float *AF = (float*) malloc(sizeof(float)*n*4*6);
	float A[6][2*n];
	Body *bodies = new Body[n];
	CudaInitialize(bs, state, n, bodies);
	

	RecDCA(bodies, n, 0, AF);

	for(int c = 1 ; c<2*n; c++)
	{
		for(int r = 0; r<6;r++)
		{
			A[r][c]=AF[2*c+r*4*n];
		}
	}	

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

		
	}
}
		
			

