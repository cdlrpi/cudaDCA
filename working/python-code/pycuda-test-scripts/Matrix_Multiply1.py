#Test script to multiply two matrices

#imports
import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule
import numpy as np
from timeit import default_timer

#Code to be sent to the GPU
#This function assumes a grid size equal to the size of the matrix, 
#where each block holds an array the length of one of the matrices sides
#Matrix multiplication consists mainly of multiplying two numbers together,
#and adding the correct products together.  This function gets all of the products
#at once, but suffers in that all of these products can not be added together simultaneously
#making the threads forced to "take turns" adding their number to the answer matrix
mod = SourceModule("""
#include <cuda.h>
__global__ void matrix_Multiply(float *A, float *B, float *C)
{
	const int i = blockIdx.x+blockIdx.y*gridDim.x; //index of the matrices related to thread position
	const int row = blockIdx.y; //The row the thread will be working with
	const int col = blockIdx.x; //The collumn the thread will be working with
	const int num = threadIdx.x; //The numbers to be multiplied together depend on the thread's x index
        const int r=row*gridDim.x+num;//The index of the desired number in the first matrix
        const int c=col+num*gridDim.x;//The index of the desired number in the second matrix
	C[i] = 0;//ensure C is initialized to zero
	
//Now the index of every number needed to complete the multiplication is found
	for(int k=0;k<blockDim.x;k++)//loop to cycle through threads
	{
		if(num == k)//The product is only added to C if it is the thread's turn
		{	    //This causes extreme slowdowns but prevents multiple threads
			    //From writing to the same location at the same time
		
			C[i]= C[i]+(A[r]*B[c]);//multiply the values at the already found indices and add them to C
		}
		__syncthreads();//Wait for the previous write to finish before attempting another write
	}
}
""")
#Create Random Matrices to Test the Function
matrix_size = 100
matrix_Multiply = mod.get_function("matrix_Multiply")
A = np.random.randn(matrix_size,matrix_size).astype(np.float32)
B = np.random.randn(matrix_size,matrix_size).astype(np.float32)
C = np.zeros_like(A)

#Start a timer and time which takes longer, the pycuda function or the numpy function
start=default_timer()
D=np.dot(A,B)
dur1=default_timer()-start
start=default_timer()

matrix_Multiply(cuda.In(A),cuda.In(B), cuda.Out(C),block=(matrix_size,1 ,1),grid=(matrix_size,matrix_size))
dur2=default_timer()-start


#Print Results
sentence1="Right Answer:"
sentence3="Time:"
sentence2="Hopefully Right Answer:"


print(sentence1)
print(D)
print(sentence3)
print(dur1)
print(sentence2)
print(C)
print(sentence3)
print(dur2)

