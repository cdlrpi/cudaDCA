
#include <stdio.h>
#include<cuda.h>
__device__ void MSsetup(float Minv[6][6],float S[6][6],int row, int col, float m, float I[],int Iindex);
__device__ void Fa(float S[6][6],float w, float r[], float m,int row, int col);
__device__ void zab(float r[],float z[6][6],float Minv[6][6],float S[6][6], int row, int col);
__device__ void zaa(float r[],float z[6][6],float Minv[6][6],float S[6][6], int row, int col);
__device__ void r02(float r[], float angle,float L);
__device__ void r01(float r[], float angle,float L);
__device__ void Mat61Mult(float A[6][6], float B[6][6], float C[6][6],int row, int col);
__device__ void Mat66Mult(float A[6][6], float B[6][6], float C[6][6],int row, int col);
__device__ void printit(float A[6][6]);
__global__ void Initialize(float state[],float m[], float l[],float I[],float Zetas[],int n)
{
    //Variable Declarations
	const int i1 = blockIdx.x;
	const int row = threadIdx.y;
	const int col = threadIdx.x;
    const int Iindex = blockIdx.x*3+row*3*n+col;
    const int i11 = col+row*26*n+i1*26;
    const int i12=i11+6;
    const int i13=i12+6; //only col=0 can write
    const int i21=i13+1;
    const int i22=i21+6;
    const int i23=i22+6; //only col = 0 can write
    __shared__ float z[6][6];
    __shared__ float Minv[6][6];
    __shared__ float S[6][6];
    __shared__ float r[3];
    __shared__ float q;
    __shared__ float w;
    //////////////////q and w//////////////////////////////////////
    //To get q and w, either this:
	__syncthreads();
    int i = 0;
	q=0;
	w=0;
    while (i <= i1)
    {
        q+=state[i];
        w+=state[i+n];
        i++;
    }
	__syncthreads();
   

    ////////////////Inverse Mass Matrix and shifter setup//////////////////////////////////
    MSsetup(Minv,S,row,col,m[i1],I,Iindex);
	//printf("here");
    
    
    ////////////////////DCM/////////////////////////////////////////
    //Dont need a dcm because it is planar, so for now to avoid memory 
    //usage, Im not using it.
    //DCM(q,C)
    
    ///////////////////z11/////////////////////////////////////////
    z[row][col]=0;
	__syncthreads();
    r01(r,q,l[i1]);
 
    zaa(r,z,Minv,S,row,col);
    
    Zetas[i11]=z[row][col];
    //////////////////z12//////////////////////////////////////////
    //S is now S10, Minv is S10*Minv and z is z11
    z[row][col]=0;
	__syncthreads();
    r02(r,q,l[i1]);
    zab(r,z,Minv,S,row,col);

    Zetas[i12]=z[row][col];
    ////////////////z13///////////////////////////////////////
    //S is S02, Minv is S10*Minv, z is z12
    r01(r,q,l[i1]);
    Fa(S,w,r,m[i1],row,col);
    Mat61Mult(Minv,S,z,row,col);
    if (col==0)
    {
        Zetas[i13]=z[row][col];
    }
    	__syncthreads();
    ////////////////z22///////////////////////////////////
    z[row][col]=0;
	__syncthreads();
    MSsetup(Minv,S,row,col,m[i1],I,Iindex);
    r02(r,q,l[i1]);
    zaa(r,z,Minv,S,row,col);
    Zetas[i22]=z[row][col];
    //////////////////z21/////////////////////////////////
    //S is now S20 Minv is S20*Minv and z is z22
    z[row][col]=0;
	__syncthreads();
    r01(r,q,l[i1]);
    zab(r,z,Minv,S,row,col);
    Zetas[i21]=z[row][col];
    ///////////////////z23/////////////////////////////////////
    //S is S01, Minv is S20*Minv and z is z21
    r02(r,q,l[i1]);
    Fa(S,w,r,m[i1],row,col);
    Mat61Mult(Minv,S,z,row,col);
    if (col==0)
    {
        Zetas[i23]=z[row][col];
    }
	__syncthreads();
}
