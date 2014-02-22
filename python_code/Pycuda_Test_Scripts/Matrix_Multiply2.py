#This script uses the gpu to multiply two matrices
#This code becomes more efficient with greater matrix sizes, but the grid and block dimensions must be
#altered so the grid dimension times the block dimension is equal to the length of one side of the matrices
#imports
import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule
import numpy as np
from timeit import default_timer
import math

#Define Block Size
BLOCK_SIZE = 100

#This is the function to run on the gpu.  It utilizes a grid pattern where each thread calculates
#one value in the answer matrix
mod = SourceModule("""
#include <cuda.h>
__global__ void matrix_Mult(float *A, float *B, float *C)
{
	const int i = threadIdx.x+threadIdx.y*blockDim.x*gridDim.x+blockIdx.x*blockDim.x+blockIdx.y*blockDim.x*blockDim.y*gridDim.x;             //This finds where in the answer matrix the thread will put it's solution
	const int row = threadIdx.y+blockIdx.y*blockDim.y;    //The row the thread is in
	const int col = threadIdx.x + blockIdx.x*blockDim.x;  //The column the thread is in

	float CC = 0;       //Create a temporary variable so global memory isn't written to each loop.
        int r=row*blockDim.x*gridDim.x;     //  The index of the beginning of the row the thread is in
	C[i]=0;          //reset the answer matrix
	
	//This loop cycles through the rows of A and the Columns of B, multiplying and adding the products
	for(int c=col;c<blockDim.x*gridDim.x*gridDim.y*blockDim.y;c=c+blockDim.x*gridDim.x)
	{
		CC+=(A[r]*B[c]); //add to the temporary variable stored locally to avoid the global write (the global reads still make this part slow)
		r=r+1;//increment to the next value in the row
	}
C[i]=CC; // Every thread saves its found value to the answer matrix.
}
"""% { 
    'BLOCK_SIZE': BLOCK_SIZE, #This is how you define values in the gpu code
    })

#get the function
matrix_Mult = mod.get_function("matrix_Mult")

#Create random matrices and compare and time the results
A = np.random.randn(BLOCK_SIZE,BLOCK_SIZE).astype(np.float32)
B = np.random.randn(BLOCK_SIZE,BLOCK_SIZE).astype(np.float32)
C = np.zeros_like(A)
start=default_timer()
D=np.dot(A,B)
dur1=default_timer()-start
start=default_timer()
matrix_Mult(cuda.In(A),cuda.In(B), cuda.Out(C),block=(10,10,1),grid=(10,10,1))
dur2=default_timer()-start
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


