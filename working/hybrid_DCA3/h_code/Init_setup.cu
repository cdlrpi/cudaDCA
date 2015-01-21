#include "../funct_bin/classes.h"
__global__ void Initialize_gpu(double state[],double m[], double l[],double I[],double Zetas[],int n);
void printa(double a[], int n);
#include <iostream>
#include <stdio.h>
void CudaInitialize(double m[], double l[], double I[], double x[], int n, double Zs[])
{
	int x_size = 2*n;
	int I_size = 3*n*3;
	int z_size = 26*6*n;
	double *d_x, *d_m, *d_l, *d_I, *d_zs;
	
	// Allocate and Load M and N to device memor

	cudaMalloc(&d_x,x_size*sizeof(double));
	cudaMalloc(&d_m, n*sizeof(double));
	cudaMalloc(&d_l, n*sizeof(double));
	cudaMalloc(&d_I, I_size*sizeof(double));
	cudaMalloc(&d_zs, z_size*sizeof(double));
	cudaMemcpy(d_x, x, x_size*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_m, m, n*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_I, I, I_size*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_l, l, n*sizeof(double), cudaMemcpyHostToDevice);
	dim3 dimBlock(6, 6,1);
	dim3 dimGrid(n,1,1);

	Initialize_gpu<<<dimGrid, dimBlock>>>(d_x, d_m, d_l, d_I, d_zs, n);
	cudaMemcpy(Zs, d_zs, z_size*sizeof(double), cudaMemcpyDeviceToHost);
	
	cudaFree(d_x);
	cudaFree(d_m);
	cudaFree(d_l);
	cudaFree(d_I);
	cudaFree(d_zs);
}
