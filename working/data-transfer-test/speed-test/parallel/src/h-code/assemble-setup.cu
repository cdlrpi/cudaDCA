
#include <iostream>
#include <stdio.h>
__global__ void Assemble_gpu(double oldZs[],double newZs[],double Xinvs[], int nn,int num);
void printm(double A[6][6]);
void printa(double A[], int x);


void cudaAssemble(double Zs[],double Xs[], int num, double nZs[], double nXs[], int odd, int newlen, float times[])
{	
	double *d_ozs,*d_nzs,*d_Xinv;
	float time1; 
	float time2;

	cudaEvent_t beginEvent1;
	cudaEvent_t endEvent1;
	cudaEvent_t beginEvent2;
	cudaEvent_t endEvent2;


	

	


	int gpulen=num-odd;
	int numBlocks = (int) (gpulen/2);
	
	cudaEventCreate( &beginEvent1 );
	cudaEventCreate( &endEvent1 );
	cudaEventRecord( beginEvent1, 0 );
	cudaMalloc(&d_ozs,sizeof(double)*(num)*6*26);
	cudaMalloc(&d_nzs,sizeof(double)*(newlen)*6*26);
	cudaMalloc(&d_Xinv,sizeof(double)*(newlen)*5*5);	
	cudaMemcpy(d_ozs, Zs, sizeof(double)*(num)*6*26, cudaMemcpyHostToDevice);
	cudaEventRecord( endEvent1, 0 );
	cudaEventSynchronize( endEvent1 );
	cudaEventElapsedTime( &time1, beginEvent1, endEvent1 );

	dim3 dimBlock(6, 6,1);
	dim3 dimGrid(numBlocks,1,1);
	cudaDeviceSynchronize();
	Assemble_gpu<<<dimGrid, dimBlock>>>(d_ozs, d_nzs, d_Xinv, newlen,num);	
	cudaDeviceSynchronize();


	cudaEventCreate( &beginEvent2 );
	cudaEventCreate( &endEvent2 );
	cudaEventRecord( beginEvent2, 0 );
	cudaMemcpy(nZs, d_nzs,sizeof(double)*(newlen)*6*26, cudaMemcpyDeviceToHost);
	cudaMemcpy(nXs, d_Xinv,sizeof(double)*(newlen)*5*5, cudaMemcpyDeviceToHost);
	cudaEventRecord( endEvent2, 0 );
	cudaEventSynchronize( endEvent2 );
	cudaEventElapsedTime( &time2, beginEvent2, endEvent2 );
	
	
	if (odd==1)
	{
		for(int r = 0; r<6; r++)
		{
			for(int c=0; c<6; c++)
			{
				nZs[c+r*newlen*26+(newlen-1)*26]=Zs[c+r*num*26+(num-1)*26];	//z11
				nZs[c+r*newlen*26+(newlen-1)*26+6]=Zs[c+r*num*26+(num-1)*26+6];	//z12
				nZs[c+r*newlen*26+(newlen-1)*26+13]=Zs[c+r*num*26+(num-1)*26+13];	//z21
				nZs[c+r*newlen*26+(newlen-1)*26+19]=Zs[c+r*num*26+(num-1)*26+19];	//z22				
				if(r != 5 && c!=5)
				{
					nXs[c+r*newlen*5+(newlen-1)*5]=Xs[c+r*num*5+(num-1)*5];
				}
			}
			nZs[r*newlen*26+(newlen-1)*26+12]=Zs[r*num*26+(num-1)*26+12];	//z13
			nZs[r*newlen*26+(newlen-1)*26+25]=Zs[r*num*26+(num-1)*26+25];	//z23
		}
	}
	//printa(nZs, newlen*6*26);
	
	times[0] += time1+time2;
	cudaFree(d_ozs);
	cudaFree(d_nzs);
	cudaFree(d_Xinv);

}

