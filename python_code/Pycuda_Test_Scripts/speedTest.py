#This script test the time it takes for various pycuda calls

#imports
import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule
import numpy as np
import time as tt
import math

#functions to make timing easier
def start():
	return tt.time()
def end(j):
	return tt.time()-j

#Functions to be run on the GPU.  These functions attempt to highlight specific gpu data tasks
kernel="""

__global__ void writetoglobal(float A) //This function writes to a global variable repeatedly
{
	for(int x=0;x<100000;x++)
	{
		A=1.4;
	}
}

__global__ void globaltoglobal(float A,float B) //Read from a global and write to a global repeatedly
{
	for(int x = 0;x<100000;x++)
	{
		B = A;
		
	}
}

__global__ void globaltoshared(float B)  //Read from a global and write to a shared repeatedly
{
	__shared__ float x;
	for(int t = 0;t<100000;t++)
	{
		x=B;
	}
}

__global__ void globaltolocal(float D)  //Read from a global and write to a local repeatedly
{
	float x;
	for(int t=0;t<100000;t++)
	{
		x=D;
	}
}

__global__ void sharedtoglobal(float C) //Read from Shared and write to global
{
	__shared__ float C_shared;
	C_shared = 5.5;
	for(int x = 0;x<100000;x++)
	{
		C = C_shared;
	}
}
__global__ void sharedtoshared(void)  //Read from shared and write to shared
{
	__shared__ float D_shared;
	__shared__ float x;
	for(int t =0;t<100000;t++)
	{
		D_shared = x;
	}
}
__global__ void sharedtolocal(void)  //Read from shared and write to local
{
	__shared__ float S;
	S=6.6;
	float x;
	for(int t =0;t<100000;t++)
	{
		x = S;
	}
}

__global__ void localtoglobal(float A) //Read from local and write to global
{	
	float x = 5.5;
	for(int t=0;t<100000;t++)
	{
		A=x;
	}
}

__global__ void localtoshared(void) //Read from local and write to shared
{
	float x=9.9;
	__shared__ float SS;
	for(int t=0;t<100000;t++)
	{
		SS=x;
	}
}

__global__ void localtolocal(void) //Read from local and write to local
{
	float x;
	float y = 1.1;
	for(int t=0;t<100000;t++)
	{
		x=y;
	}
}
__global__ void latency(void) //A dummy function that requires no input
{
}	
__global__ void memDummy1(float* F) //A dummy function that require an array of floats
{
}
__global__ void memDummy2(float F) //A dummy function that requires a single float
{
}
__global__ void memDummy3(int x) // A dummy function that requires an integer
{
}
"""

#Compile the code and get the functions
mod = SourceModule(kernel)
writetoglobal=mod.get_function("writetoglobal")
globaltoglobal= mod.get_function("globaltoglobal")
globaltoshared=mod.get_function("globaltoshared")
globaltolocal=mod.get_function("globaltolocal")
sharedtoglobal=mod.get_function("sharedtoglobal")
sharedtoshared=mod.get_function("sharedtoshared")
sharedtolocal=mod.get_function("sharedtolocal")
localtoglobal=mod.get_function("localtoglobal")
localtoshared=mod.get_function("localtoshared")
localtolocal=mod.get_function("localtolocal")
memDummy1=mod.get_function("memDummy1")
memDummy2=mod.get_function("memDummy2")
memDummy3=mod.get_function("memDummy3")
latency = mod.get_function("latency")

#Create a few floats to send to the gpu to test
A=np.array([6.7],dtype=np.float32)
B=np.float32(7.7)
C=np.float32(8.7)
D=np.float32(9.7)
E=np.float32(17)
F=np.float32(8.4)
G=np.float32(7.9)
H=np.float32(5.4)
arr=np.random.random_sample((100000,1))
natlat1=0.0

#The following code calls each of the functions above that involve reading and writing from specific memory locations
#All of the function calls are timed
#write to global
i=start()
writetoglobal(cuda.In(B),block=(1024,1,1),grid=(1,1,1))
writetoglobaldur=end(i)
#global to global
i=start()
globaltoglobal(cuda.In(C),block=(1024,1,1),grid=(1,1,1))
globaltoglobaldur=end(i)

#global to shared
i=start()
globaltoshared(cuda.In(D),block=(1024,1,1),grid=(1,1,1))
globaltoshareddur=end(i)

#global to local
i=start()
globaltolocal(cuda.In(E),block=(1024,1,1),grid=(1,1,1))
globaltolocaldur=end(i)

#shared to global
i=start()
sharedtoglobal(cuda.In(F),block=(1024,1,1),grid=(1,1,1))
sharedtoglobaldur=end(i)

#shared to shared
i=start()
sharedtoshared(block=(1024,1,1),grid=(1,1,1))
sharedtoshareddur=end(i)
#shared to local
i=start()
sharedtolocal(block=(1024,1,1),grid=(1,1,1))
sharedtolocaldur=end(i)
#local to global
i=start()
localtoglobal(cuda.In(G),block=(1024,1,1),grid=(1,1,1))
localtoglobaldur=end(i)

#local to shared
i=start()
localtoshared(block=(1024,1,1),grid=(1,1,1))
localtoshareddur=end(i)

#local to local
i=start()
localtolocal(block=(1024,1,1),grid=(1,1,1))
localtolocaldur=end(i)

#get the latency due to calling the function when no data is transfered
i=0.0
for x in range(0,100):
	i=start() 
	latency(block=(1,1,1),grid=(1,1,1))#with small block size
	i=end(i)
	natlat1+=i
natlat1=natlat1/100 #save the average time the call takes

i=0.0
natlat2=0;
for x in range(0,100):
	i=start()
	latency(block=(1024,1,1),grid=(1,1,1))#large block size
	i=end(i)
	natlat2+=i
natlat2=natlat2/100 #save the average time the call takes


#latency from data send (only 1 float)
i=0.0
sendlat=0.0;
for x in range(0,100):
	i=start()
	memDummy2(cuda.In(H),block=(1,1,1),grid=(1,1,1))
	i=end(i)
	
	sendlat+=i
sendlat/=100 #Save the average time it takes to send a single float

#memory test-1000000 floats,1 array
onefloat=0.0
for x in range(0,100):
	i=start()
	memDummy1(cuda.In(arr),block=(1,1,1),grid=(1,1,1))
	onefloat+=(end(i)-sendlat)/100000
onefloat/=100 #attempt to find the transfer time per float in the array by subtracting the time it takes 
		#to send a single float, hopefully eliminating the latency of sending any data.
sendlat=sendlat-onefloat #Adjust the time it takes to send data by subtracting the time it took for the single
				#float to be send



#print results
print()
print("Latency of Function Call (1 thread): ",natlat1)
print("Latency of Function Call (1024 threads): ",natlat2)
print("Latency of Function Call with Memory Send: ",sendlat)
print("Time per float transfered in an array: ",onefloat)

print()
print("Read From | Write To | Time")
print("  Global  |  Global  |",(globaltoglobaldur-(sendlat)-(2*onefloat))*10000)
print("  Global  |  Shared  |",10000*(globaltoshareddur-sendlat-onefloat))
print("  Global  |  Local   |",10000*(globaltolocaldur-sendlat-onefloat))
print("  Shared  |  Global  |",10000*(sharedtoglobaldur-sendlat-onefloat))
print("  Shared  |  Shared  |",10000*(sharedtoshareddur-natlat2))
print("  Shared  |  Local   |",10000*(sharedtolocaldur-natlat2))
print("  Local   |  Global  |",10000*(localtoglobaldur-sendlat-onefloat))
print("  Local   |  Shared  |",10000*(localtoshareddur-natlat2))
print("  Local   |  Local   |",10000*(localtolocaldur-natlat2))
print()
print("Global write time (no read): ",(writetoglobaldur-onefloat-sendlat)*10000)









