#include <iostream>
__global__ void Disassemble_gpu(double Xinv[],double Zs[],double oldAF[], double newAF[], int numBlocks, int lesslen);
void printm(double A[6][6]);
void printa(double A[], int x);


void cudaDisassemble(double OldAF[], double Zs[], double Xs[],double nZs[], double nXs[], int odd, int morelen, int lesslen, double AF[], float times[])
{	
	double *d_oldAF, *d_newAF,*d_Xinv;
	double *d_zs;
	float time1; 
	float time2;

	cudaEvent_t beginEvent1;
	cudaEvent_t endEvent1;
	cudaEvent_t beginEvent2;
	cudaEvent_t endEvent2;
	
	int newlen = (int) morelen*4;
	int numBlocks = (int) morelen-lesslen;
 
 	cudaEventCreate( &beginEvent1 );
	cudaEventCreate( &endEvent1 );
	cudaEventRecord( beginEvent1, 0 );
	cudaMalloc(&d_zs,sizeof(double)*(morelen)*6*26);
	cudaMalloc(&d_newAF,sizeof(double)*(newlen)*6);
	cudaMalloc(&d_oldAF,sizeof(double)*(lesslen)*4*6);
	cudaMalloc(&d_Xinv,sizeof(double)*(lesslen)*5*5);	
	cudaMemcpy(d_zs, Zs, sizeof(double)*(morelen)*6*26, cudaMemcpyHostToDevice);
	cudaMemcpy(d_Xinv, nXs, sizeof(double)*(lesslen)*5*5, cudaMemcpyHostToDevice);
	cudaMemcpy(d_oldAF, OldAF, sizeof(double)*(lesslen)*4*6, cudaMemcpyHostToDevice);
	cudaEventRecord( endEvent1, 0 );
	cudaEventSynchronize( endEvent1 );
	cudaEventElapsedTime( &time1, beginEvent1, endEvent1 );

	dim3 dimBlock(6, 6,1);
	dim3 dimGrid(numBlocks,1,1);
	
	Disassemble_gpu<<<dimGrid, dimBlock>>>(d_Xinv, d_zs, d_oldAF, d_newAF, morelen,lesslen);
	cudaEventCreate( &beginEvent2 );
	cudaEventCreate( &endEvent2 );
	cudaEventRecord( beginEvent2, 0 );
	cudaMemcpy(AF, d_newAF,sizeof(double)*(newlen)*6, cudaMemcpyDeviceToHost);
	cudaEventRecord( endEvent2, 0 );
	cudaEventSynchronize( endEvent2 );
	cudaEventElapsedTime( &time2, beginEvent2, endEvent2 );
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
	times[0] += time1+time2;
	cudaFree(d_zs);
	cudaFree(d_newAF);
	cudaFree(d_oldAF);
	cudaFree(d_Xinv);
	
	

}

	

