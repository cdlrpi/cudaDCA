//
//Cuda DCA Driver
//
//So the next thing to do is to load m_gpu, l_gpu, and I_gpu onto the gpu.  
//Then I write something to get initialize ready and interperate the output
//and I check if it matches the initialization pycuda function.
//Then I do the same for Assembly and see if it matches the final assembly
//values from pycuda.  Then we try disassembly and compare again.  Then I 
//get the final acceleration values and check that those are correct.  And
//then I can finally try out the integrator and hope the whole thing works!
#include <malloc.h>
#include "h_code/MainHeader.h"
#include <cuda.h>
#include <cuda_runtime.h>
//Main function
int main()
{
	//Variable Declarations
	InitBody* bodies;
	Joint* joints;
	float *inits;
	float *Time;
	float *Y;
	float **X;
	
	//System Setup
	int n=5;
	bodies = new InitBody[n];
	joints = new Joint[n];
	inits = new float[2*n];

	//Time Setup
	float tstep=.01;
	float tfinal = 4.0;
	int tlen = (int) floor(tfinal/tstep);
	Time = new float[tlen];

	//Allocation of memory for the matrix holding the solution 
	Y = new float[2*n];
	X = new float*[tlen];
	for ( int i = 0; i < tlen ;i++)
	{
		X[i]= new float[2*n];
	}
	//Initialize lengths and masses to 1
	pend_init(bodies,n,1.0,1.0);
	
	//Initialize the position and velocity to match dropping an n link pendulum 
	//where all links are parallel to the horizontal and velocity is 0
	full_drop(inits,n);

	//Initialize the time array
	Time_init(Time,tstep,tfinal);


	
	//practice call to the integrator
	//RK_45(inits, tstep, n , bodies, joints, Y);

	DCAhelp(inits,bodies, joints, n, Y);
	
	delete[] inits;
	delete[] Time;
	delete[] Y;
	delete[] X;

	return 0;

}


