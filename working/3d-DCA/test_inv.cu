#include <iostream>
#include <stdio.h>
__global__ void Invert_Matrix(double mat[],double eye[]);
int main(void)
{
	double *d_A, *d_B;
	double mat3[36];
	double mat[6][6]={{0.060471,0.291984,0.372410,0.052677,0.417744,0.698106},{0.399258,0.431651,0.198118,0.737858,0.983052,0.666528},{0.526876,0.015487,0.489688,0.269119,0.301455,0.178132},{0.416799,0.984064,0.339493,0.422836,0.701099,0.128014},{0.656860,0.167168,0.951630,0.547871,0.666339,0.999080},{0.627973,0.106216,0.920332,0.942737,0.539126,0.171121}};
	double mat2[36];
	for(int i =0; i<6; i++)
	{
		for(int j=0; j<6; j++)
		{
			mat2[i*6+j]=mat[i][j];
		}
	}
	double eye[6][6];
	double eye2[36];
	for(int i =0; i<6; i++)
	{
		for(int j=0; j<6; j++)
		{
			if(j==i)
			{
				eye2[i*6+j]=1;
				eye[i][i]=1;
			}	
			else
			{
				eye2[i*6+j]=0;
				eye[i][j]=0;
			}
		}
	}
	double temp1;
	double temp2;
	std::cout<<"haha"<<std::endl;
	double ar[6*6] = {5,2,3,4,5,6,7,8,1,10,11,12,13,14,15,16,5,18,19,20,21,3,23,24,25,26,27,28,29,30,31,32,33,34,35,36};
	std::cout<<"blah"<<std::endl;
/*
	for(int r =0; r<6; r++)
	{
		for(int c=0; c<6; c++)
		{
			mat[r][c]= ar[r*6+c];
		}
		std::cout<<"hey"<<std::endl;
	}
	std::cout<<"what the?"<<std::endl;
*/
for(int x =0; x<6; x++)
		{
			for(int y =0; y<6; y++)
			{
				std::cout<<mat[x][y]<<'\t';
			}
		std::cout<<std::endl;
		}
	for(int r =0; r<6; r++)
	{
		for(int x =0; x<6; x++)
		{
			for(int y =0; y<6; y++)
			{
				std::cout<<eye[x][y]<<'\t';
			}
		std::cout<<std::endl;
		}
		std::cout<<std::endl;
		temp1=mat[r][r];
		for(int c =0; c<6; c++)
		{
			mat[r][c]=mat[r][c]/temp1;
			eye[r][c]=eye[r][c]/temp1;
			//std::cout<<'x'<<temp1<<std::endl;
			//std::cout<<'y'<<r<<'c'<<c<<std::endl;
			//std::cout<<mat[r][c]<<std::endl;
		}
		for(int i=0; i<6; i++)
		{
			if(i!=r)
			{
				temp2=mat[i][r];
				for(int c =0; c<6; c++)
				{
					mat[i][c]-=(mat[r][c]*temp2);
					eye[i][c]-=(eye[r][c]*temp2);
				}
			}
		}
	}
	for(int r =0; r<6; r++)
	{
		for(int c =0; c<6; c++)
		{
			std::cout<<mat[r][c]<<'\t';
		}
		std::cout<<std::endl;
	}
	for(int r =0; r<6; r++)
	{
		for(int c =0; c<6; c++)
		{
			std::cout<<eye[r][c]<<'\t';
		}
		std::cout<<std::endl;
	}
	cudaMalloc(&d_A,sizeof(double)*36);	
	cudaMalloc(&d_B,sizeof(double)*36);
	cudaMemcpy(d_A,mat2,sizeof(double)*36,cudaMemcpyHostToDevice);
	cudaMemcpy(d_B,eye2,sizeof(double)*36,cudaMemcpyHostToDevice);
	dim3 dimBlock(6, 6,1);
	dim3 dimGrid(1,1,1);
	cudaDeviceSynchronize();
	Invert_Matrix<<<dimGrid, dimBlock>>>(d_A,d_B);	
	cudaDeviceSynchronize();
	cudaMemcpy(mat3, d_A,sizeof(double)*36, cudaMemcpyDeviceToHost);
	cudaMemcpy(eye2, d_B,sizeof(double)*36, cudaMemcpyDeviceToHost);
	std::cout<<std::endl<<std::endl<<std::endl;
	for(int r =0; r<6; r++)
	{
		for(int c =0; c<6; c++)
		{
			std::cout<<mat3[r*6+c]<<'\t';
		}
		std::cout<<std::endl;
	}
}

__global__ void Invert_Matrix(double A[],double B[])
{	
	const int row = threadIdx.y;
	const int col= threadIdx.x;
	__shared__ double temp1;
	__shared__ double temp2;
	__shared__ double mat[6][6];
	__shared__ double eye[6][6];
	mat[row][col]=A[row*6+col];
	eye[row][col]=B[row*6+col];
	__syncthreads();
	if(row == col)
	{
		for(int r =0; r<6; r++)
		{
			if(col==0)
			{
				
				temp1=mat[r][r];
			}
			__syncthreads();
			mat[r][col]/=temp1;
			eye[r][col]/=temp1;

			
			for(int i=0; i<6; i++)
			{
				if(i!=r)
				{
					if(col==0)
					{
						temp2=mat[i][r];
					}
					__syncthreads();
					mat[i][col]-=(mat[r][col]*temp2);
					eye[i][col]-=(eye[r][col]*temp2);
				}
			}
		}
	}
	
	A[row*6+col]=mat[row][col];
	B[row*6+col]=eye[row][col];	
	__syncthreads();
	printf("xx%fxx",A[row*6+col]);
	if(row==0 && col==0)
	{
		for(int j =0; j<6; j++)
				{
					for(int i =0; i<6; i++)
					{
						printf("%f\t",A[j*6+i]);
					}
					printf("\n\n\n");
				}
	}
}			
	
