#include "classes.h"
#include <iostream>
__global__ void Disassemble(float Xinv[],float Zs[],float oldAF[], float newAF[], int numBlocks);
void printm(float A[6][6]);
void printa(float A[], int x);
void cudaDisassemble(float OldAF[], Body *morebds, Body *lessbds, int odd, int morelen, int lesslen, float AF[])
{	
	float *d_oldAF, *d_newAF,*d_Xinv, *zs, *newAF, *Xinv, *oldAF;
	float *d_zs;
	int r11,r12,r13,r21,r22,r23,rXinv;
	int newlen = (int) morelen*4;
	int gpu_len = (int) (newlen-(4*odd));
	int numBlocks = (int) morelen-lesslen;
	zs = (float*) malloc(sizeof(float)*(numBlocks*2)*6*26);
	newAF = (float*) malloc(sizeof(float)*(gpu_len)*6);
	oldAF = (float*) malloc(sizeof(float)*(lesslen-odd)*4*6);
	Xinv = (float*) malloc(sizeof(float)*(numBlocks)*5*5);

	
	for(int c = 0; c<(lesslen-odd)*4; c++)
	{
		for(int r = 0; r<6; r++)
		{
			oldAF[c+r*(lesslen-odd)*4]= OldAF[c+r*lesslen*4];
		}
	}
	for(int i = 0; i<numBlocks*2; i++)
	{
		for(int c = 0; c<6; c++)
		{
			for(int r = 0; r<6; r++)
			{	r11=26*(numBlocks*2)*r+c+i*26;
				r12=r11+6;
				r13=r12+6;
				r21=r13+1;
				r22=r21+6;
				r23=r22+6;
				
				zs[r11] = morebds[i].z11[r][c];
				zs[r12] = morebds[i].z12[r][c];
				zs[r21] = morebds[i].z21[r][c];
				zs[r22] = morebds[i].z22[r][c];
				
				if(c==0)
				{
					zs[r13] = morebds[i].z13[r];
					zs[r23] = morebds[i].z23[r];
				}
			}
		}
	}

	for(int i = 0; i<numBlocks; i++)
	{
		for(int c = 0; c<5; c++)
		{
			for(int r = 0; r<5; r++)
			{	
				rXinv = i*5+c+r*5*numBlocks;
				Xinv[rXinv]= lessbds[i].Xinv[r][c];
			}
		}
	}

	
   
	std::cout<<cudaMalloc(&d_zs,sizeof(float)*(numBlocks*2)*6*26);
	std::cout<<cudaDeviceSynchronize();
	std::cout<<cudaMalloc(&d_newAF,sizeof(float)*(gpu_len)*6);
	std::cout<<cudaDeviceSynchronize();
	std::cout<<cudaMalloc(&d_oldAF,sizeof(float)*(lesslen-odd)*4*6);
	std::cout<<cudaDeviceSynchronize();
	std::cout<<cudaMalloc(&d_Xinv,sizeof(float)*(numBlocks)*5*5);	
std::cout<<cudaDeviceSynchronize();
	std::cout<<cudaMemcpy(d_zs, zs, sizeof(float)*(numBlocks*2)*6*26, cudaMemcpyHostToDevice);
std::cout<<cudaDeviceSynchronize();
	std::cout<<cudaMemcpy(d_Xinv, Xinv, sizeof(float)*(numBlocks)*5*5, cudaMemcpyHostToDevice);
std::cout<<cudaDeviceSynchronize();
	std::cout<<cudaMemcpy(d_oldAF, oldAF, sizeof(float)*(lesslen-odd)*4*6, cudaMemcpyHostToDevice);

	dim3 dimBlock(6, 6,1);
	dim3 dimGrid(numBlocks,1,1);
	std::cout<<cudaDeviceSynchronize();
	Disassemble<<<dimGrid, dimBlock>>>(d_Xinv, d_zs, d_oldAF, d_newAF, numBlocks);
	std::cout<<cudaDeviceSynchronize();
	std::cout<<cudaMemcpy(newAF, d_newAF,sizeof(float)*(gpu_len)*6, cudaMemcpyDeviceToHost);
	std::cout<<cudaDeviceSynchronize();
	for(int c = 0; c<(morelen)*4; c++)
	{
		for (int r = 0; r<6;r++)
		{	
			
			AF[r*morelen*4+c] = newAF[c+r*(gpu_len)];
		
		}
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
	
	std::cout<<cudaFree(d_zs);
	std::cout<<cudaFree(d_newAF);
	std::cout<<cudaFree(d_oldAF);
	std::cout<<cudaFree(d_Xinv);
	
	free(zs);
	free(newAF);
	free(oldAF);
	free(Xinv);

}

	
