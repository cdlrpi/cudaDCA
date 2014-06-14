#include "classes.h"
void matmult61(float A[6][6], float B[6], float C[6]);



void solve_BCs(Body *bodies, float AF[])
{
	float* Fc1 = (float*) malloc(sizeof(float)*6);
	float* Fc2 = (float*) malloc(sizeof(float)*6);
	float* A1 = (float*) malloc(sizeof(float)*6);
	float* A2 = (float*) malloc(sizeof(float)*6);
	float* temp = (float*)malloc(sizeof(float)*6);
	float M[3][4];
	
	for(int c = 0; c<3; c++)
	{
		for(int r = 0; r<3; r++)
		{
			M[r][c] = bodies[0].z11[r+3][c+3];
			if (c==0)
			{
				M[r][3] = -1*bodies[0].z13[r+3];
			}
		}
	}
	for(int s=0; s<3;s++)
	{
		for(int j=0; j<4;j++)
		{
			M[s][j]=M[s][j]/M[s][s];
		}
		for(int k =0; k<3; k++)
		{
			if( s!=k)
			{
				for(int l=0; l<4; l++)
				{
					M[k][l]-=(M[s][l]*M[k][s]);
				}
			}
		}
	}
	for(int s=0; s<6;s++)
	{
		Fc2[s] = 0;
		if (s>=3)
		{
			Fc1[s] = M[s-3][3];
		}
		else
		{
			Fc1[s] = 0;
		}
	}
	matmult61(bodies[0].z11, Fc1, temp);
	for(int s = 0; s<6; s++)
	{
		A1[s]=temp[s]+bodies[0].z13[s];
	}
	matmult61(bodies[0].z21, Fc1, temp);
	for(int s = 0; s<6; s++)
	{
		A2[s]=temp[s]+bodies[0].z23[s];
	}
	
	
	for(int r =0; r<6; r++)
	{
		AF[r*4]=A1[r];
		AF[r*4+1]=Fc1[r];
		AF[r*4+2]=A2[r];
		AF[r*4+3]=Fc2[r];
	}
}

	

	

