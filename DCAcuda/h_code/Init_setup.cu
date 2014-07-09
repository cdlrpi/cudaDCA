#include "classes.h"
__global__ void Initialize(double state[],double m[], double l[],double I[],double Zetas[],int n);
void printa(double a[], int n);
#include <iostream>
#include <stdio.h>
void CudaInitialize(InitBody* oldbds, double x[], int n, Body *newbds)
{
	int x_size = 2*n;
	int I_size = 3*n*3;
	int z_size = 26*6*n;
	int r11,r12,r13,r21,r22,r23;
	double *x_gpu= (double*)malloc(x_size*sizeof(double));
	double *m_gpu =(double*)malloc(n*sizeof(double));
	double *l_gpu = (double*)malloc(n*sizeof(double));
	double *I_gpu = (double*)malloc(I_size*sizeof(double));
	double *zs_gpu = (double*)malloc(z_size*sizeof(double));
	double *d_x, *d_m, *d_l, *d_I, *d_zs;
	
	for(int i = 0; i<n; i++)
	{
		x_gpu[i] = x[i];
		x_gpu[i+n] = x[i + n];
		m_gpu[i] = oldbds[i].m;
		l_gpu[i]= oldbds[i].l;
		for(int c = 0; c<3; c++)
		{
			for(int r = 0; r<3; r++)
			{
				I_gpu[3*n*r+c+i*3]=oldbds[i].I[r][c];
			}
		}
	}

	// Allocate and Load M and N to device memor

	cudaMalloc(&d_x,x_size*sizeof(double));
	cudaMalloc(&d_m, n*sizeof(double));
	cudaMalloc(&d_l, n*sizeof(double));
	cudaMalloc(&d_I, I_size*sizeof(double));
	cudaMalloc(&d_zs, z_size*sizeof(double));
	cudaMemcpy(d_x, x_gpu, x_size*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_m, m_gpu, n*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_I, I_gpu, I_size*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_l, l_gpu, n*sizeof(double), cudaMemcpyHostToDevice);
	dim3 dimBlock(6, 6,1);
	dim3 dimGrid(n,1,1);

	Initialize<<<dimGrid, dimBlock>>>(d_x, d_m, d_l, d_I, d_zs, n);
	cudaMemcpy(zs_gpu, d_zs, z_size*sizeof(double), cudaMemcpyDeviceToHost);
	for(int i = 0; i<n; i++)
	{
		for(int c = 0; c<6; c++)
		{
			for(int r = 0; r<6; r++)
			{	r11=26*n*r+c+i*26;
				r12=r11+6;
				r13=r12+6;
				r21=r13+1;
				r22=r21+6;
				r23=r22+6;
				newbds[i].z11[r][c]=zs_gpu[r11];
				newbds[i].z12[r][c]=zs_gpu[r12];
				newbds[i].z21[r][c]=zs_gpu[r21];
				newbds[i].z22[r][c]=zs_gpu[r22];
				if(c==0)
				{
					newbds[i].z13[r]=zs_gpu[r13];
					newbds[i].z23[r]=zs_gpu[r23];
				}
			}
			
		}

	
	}


	cudaFree(d_x);
	cudaFree(d_m);
	cudaFree(d_l);
	cudaFree(d_I);
	cudaFree(d_zs);
	
free(x_gpu);
	free(m_gpu);
	free(l_gpu);
	free(I_gpu);
	free(zs_gpu);

}
