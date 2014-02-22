#This script continually sends arrays of increasing size to the gpu
#Timing each send, and graphing the time it takes to transfer vs the length of the arrays being transfered

#Imports
import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule
import numpy as np
from timeit import default_timer
import math
import matplotlib.pyplot as plt

#Code that goes to the gpu.  This function does nothing but accept input
kernel="""

__global__ void memDummy1(float* F,float* D, float* C)
{
}

"""
#compile code and get function
mod = SourceModule(kernel)
memDummy=mod.get_function("memDummy1")

#Initialize necessary variables and begin loop, collecting data
y=np.zeros((100,1))
arraylength=1
i=0
while arraylength <=1000000:
	A=B=C=np.random.random_sample((arraylength,1)).astype(np.float32)#Get three random arrays
	start=default_timer() #start the timer
	memDummy(cuda.InOut(A),cuda.In(B),cuda.In(C),block=(1,1,1),grid=(1,1,1)) #call the function
	y[i]=default_timer()-start #end the timer and save the value
	arraylength+=10000 #increment the array length
	i+=1
#plot the graph and label the axes
plt.plot(np.linspace(0,1000000,100),y)
plt.ylabel('Time to transfer')	
plt.xlabel('Length of Arrays')
plt.show()	

