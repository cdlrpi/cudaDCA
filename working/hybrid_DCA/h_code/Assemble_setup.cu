#include "../funct_bin/classes.h"
#include <iostream>
__global__ void Assemble_gpu(double oldZs[],double newZs[],double Xinvs[], int nn);
void printm(double A[6][6]);
void printa(double A[], int x);
void cudaAssemble(Body *bodies, int num, Body *newbds, int odd, int newlen)
{	
	double *d_ozs,*d_nzs,*d_Xinv, *ozs, *nzs, *Xinv;
	int r11,r12,r13,r21,r22,r23,rXinv;
	
	int gpulen=num-odd;
	int numBlocks = (int) (gpulen/2);
	ozs = (double*) malloc(sizeof(double)*(num-odd)*6*26);
	nzs = (double*) malloc(sizeof(double)*(newlen-odd)*6*26);
	Xinv = (double*) malloc(sizeof(double)*(newlen-odd)*5*5);
	
	for(int i = 0; i<gpulen; i++)
	{
		for(int c = 0; c<6; c++)
		{
			for(int r = 0; r<6; r++)
			{	r11=26*gpulen*r+c+i*26;
				r12=r11+6;
				r13=r12+6;
				r21=r13+1;
				r22=r21+6;
				r23=r22+6;
				ozs[r11] = bodies[i].z11[r][c];
				ozs[r12] = bodies[i].z12[r][c];
				ozs[r21] = bodies[i].z21[r][c];
				ozs[r22] = bodies[i].z22[r][c];

				if(c==0)
				{
					ozs[r13] = bodies[i].z13[r];
					ozs[r23] = bodies[i].z23[r];
				}
			}
		}
	}
	
	cudaMalloc(&d_ozs,sizeof(double)*(num-odd)*6*26);
	cudaMalloc(&d_nzs,sizeof(double)*(newlen-odd)*6*26);
	cudaMalloc(&d_Xinv,sizeof(double)*(newlen-odd)*5*5);	
	cudaMemcpy(d_ozs, ozs, sizeof(double)*(num-odd)*6*26, cudaMemcpyHostToDevice);

	dim3 dimBlock(6, 6,1);
	dim3 dimGrid(numBlocks,1,1);
	Assemble_gpu<<<dimGrid, dimBlock>>>(d_ozs, d_nzs, d_Xinv, numBlocks);	
	
	if (odd == 1)
	{
		for(int r = 0; r<6; r++)
		{
			for(int c=0; c<6; c++)
			{
				newbds[newlen-1].z11[r][c]=bodies[num-1].z11[r][c];
				newbds[newlen-1].z12[r][c]=bodies[num-1].z12[r][c];
				newbds[newlen-1].z21[r][c]=bodies[num-1].z21[r][c];
				newbds[newlen-1].z22[r][c]=bodies[num-1].z22[r][c];
				
				if(r != 5 && c!=5)
				{
					newbds[newlen-1].Xinv[r][c]=bodies[num-1].Xinv[r][c];
				}
			}
			newbds[newlen-1].z13[r]=bodies[num-1].z13[r];
			newbds[newlen-1].z23[r]=bodies[num-1].z23[r];
		}
	}

	cudaMemcpy(nzs, d_nzs,sizeof(double)*(newlen-odd)*6*26, cudaMemcpyDeviceToHost);
	cudaMemcpy(Xinv, d_Xinv,sizeof(double)*(newlen-odd)*5*5, cudaMemcpyDeviceToHost);

	for(int i = 0; i<newlen; i++)
	{
		for(int c = 0; c<6; c++)
		{
			for(int r = 0; r<6; r++)
			{	
				if(!(odd ==1 && i==newlen-1))
				{
					r11=26*(newlen-odd)*r+c+i*26;
					r12=r11+6;
					r13=r12+6;
					r21=r13+1;
					r22=r21+6;
					r23=r22+6;
					rXinv = 5*i+c+r*(newlen-odd)*5;
					newbds[i].z11[r][c] = nzs[r11];
					newbds[i].z12[r][c] = nzs[r12]; 
					newbds[i].z21[r][c] = nzs[r21]; 
					newbds[i].z22[r][c] = nzs[r22]; 
					
					if(r!=5 && c!=5)
					{
						newbds[i].Xinv[r][c]= Xinv[rXinv];
					}
					if(c==0)
					{	
						newbds[i].z13[r] = nzs[r13]; 
						newbds[i].z23[r] = nzs[r23]; 
					}
				}
			}
		}
	
	}

	cudaFree(d_ozs);
	cudaFree(d_nzs);
	cudaFree(d_Xinv);
	
free(ozs);
	free(nzs);
	free(Xinv);

}
