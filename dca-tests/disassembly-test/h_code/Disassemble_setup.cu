#include "../funct_bin/classes.h"
#include <iostream>
__global__ void Disassemble_gpu(double Xinv[],double Zs[],double oldAF[], double newAF[], int numBlocks, int lesslen);
void printm(double A[6][6]);
void printa(double A[], int x);


void cudaDisassemble(double OldAF[], double Zs[], double Xs[],double nZs[], double nXs[], int odd, int morelen, int lesslen, double AF[],int data)
{	
	double *d_oldAF, *d_newAF,*d_Xinv;
	double *d_zs;
	
	int newlen = (int) morelen*4;
	int numBlocks = (int) morelen-lesslen;
 
	cudaMalloc(&d_zs,sizeof(double)*(morelen)*6*26);
	cudaMalloc(&d_newAF,sizeof(double)*(newlen)*6);
	cudaMalloc(&d_oldAF,sizeof(double)*(lesslen)*4*6);
	cudaMalloc(&d_Xinv,sizeof(double)*(lesslen)*5*5);	
	if(data)
	{
	cudaMemcpy(d_zs, Zs, sizeof(double)*(morelen)*6*26, cudaMemcpyHostToDevice);
	cudaMemcpy(d_Xinv, nXs, sizeof(double)*(lesslen)*5*5, cudaMemcpyHostToDevice);
	cudaMemcpy(d_oldAF, OldAF, sizeof(double)*(lesslen)*4*6, cudaMemcpyHostToDevice);
	}
	dim3 dimBlock(6, 6,1);
	dim3 dimGrid(numBlocks,1,1);
	
	Disassemble_gpu<<<dimGrid, dimBlock>>>(d_Xinv, d_zs, d_oldAF, d_newAF, morelen,lesslen);
	if(data)
	{
	cudaMemcpy(AF, d_newAF,sizeof(double)*(newlen)*6, cudaMemcpyDeviceToHost);
	}
	if(odd ==1)
	{
		for (int r = 0; r<6;r++)
		{
			AF[r*morelen*4+morelen*4-4]=OldAF[r*lesslen*4+lesslen*4-4];
			AF[r*morelen*4+morelen*4-3]=OldAF[r*lesslen*4+lesslen*4-3];
			AF[r*morelen*4+morelen*4-2]=OldAF[r*lesslen*4+lesslen*4-2];
			AF[r*morelen*4+morelen*4-1]=OldAF[r*lesslen*4+lesslen*4-1];
		}
	}
	
	//std::cin.get();
	
	cudaFree(d_zs);
	cudaFree(d_newAF);
	cudaFree(d_oldAF);
	cudaFree(d_Xinv);
	
	

}

	

	
