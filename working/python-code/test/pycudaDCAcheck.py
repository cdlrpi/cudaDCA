
# coding: utf-8

# In[93]:

# Pycuda DCA
#---------------
# This is a script to solve for the motion of a 
# massive rod double pendulum using DCA and pycuda

import numpy as np
import math
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import time

import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule

import MBstructs as MB
import MultiBodyFuncts as MBF




# In[94]:

# System Definition
#--------------------

#number of bodies
n=21

#Create a list of Bodies and Joints
bs=[]
js=[]

for i in range(0,n):
    bs.append(MB.Body())
    js.append(MB.Joint())

#Joint Properties
for joint in js:
    joint.PDmatrix(np.array((0,0,1,0,0,0)),1)

#Body Properties
for i in range (0,n):
    bs[i].m=1.0
    bs[i].l=1.0
    bs[i].Inertia(bs[i].m*(bs[i].l**2)/12,1,bs[i].m*(bs[i].l**2)/12)


# In[95]:

# Initial Conditions
#--------------------
# q1 and q2 are measured as the angle the rod makes with -n2 axis
q_init=np.zeros((n))
q_init[0]=np.pi/2

qdot_init=np.zeros((n))

initialconditions=np.zeros((2*n))
initialconditions[:n]=q_init
initialconditions[n:]=qdot_init

#Length of time of the simulation

Time=np.arange(0,1,.01)


# In[96]:

mod = SourceModule("""
/*This is a cuda kernal which calculates the initial zeta values for all bodies
in a pendulum chain*/
#include <cuda.h>
#include <stdio.h>

//Function Prototypes
__device__ void set_Dt_pend(float D[6][6],int row,int col);
__device__ void make_W_pend(float Xinv[6][6], float D[6][6],int row , int col);                    
__device__ void get_X_pend(float z1[6][6], float z2 [6][6], float D[6][6],int row , int col);
__device__ void set_D_pend(float D[6][6],int row,int col);
__device__ void get_W(float z1[6][6], float z2[6][6], float W[6][6],int row, int col);
__device__ void invert_X(float X[6][6], float A[6][6],int row, int col);
__device__ void PrintMat61(float A[6][6], int nn );
__device__ void PrintMat66(float A[6][6], int nn);
__device__ void MSsetup(float Minv[6][6],float S[6][6],int row, int col, float m, float I[],int Iindex,int i1);
__device__ void Fa(float S[6][6],float w, float r[], float m,int row, int col);
__device__ void transpose(float S[6][6],int row, int col);
__device__ void makeS(float S[6][6],float *r);
__device__ void zab(float r[],float z[6][6],float Minv[6][6],float S[6][6], int row, int col);
__device__ void zaa(float r[],float z[6][6],float Minv[6][6],float S[6][6], int row, int col);
__device__ void r02(float r[], float angle,float L);
__device__ void r01(float r[], float angle,float L);
__device__ void DCM(float angle, float **B);
__device__ void Mat61Mult(float A[6][6], float B[6][6], float C[6][6],int row, int col);
__device__ void Mat66Mult(float A[6][6], float B[6][6], float C[6][6],int row, int col);

__global__ void Disassemble(float Xinv[] ,float zs[], float oldAs[] ,float newAs[], int numBlocks)
{
    const int i1 = blockIdx.x;
    const int row = threadIdx.y;
    const int col = threadIdx.x;
//index to read Xinv[]
    const int rXinv=col+row*numBlocks*5+i1*5;
//indices to read zeta values
    const int r111=col+row*numBlocks*52+i1*52;
    const int r112=r111+6;
    const int r113=r112+6; //only col = 0 can read
    const int r121=r113+1;
    const int r122=r121+6;
    const int r123=r122+6; //only col= 0 can read
    
    const int r211=r123+1;
    const int r212=r211+6;
    const int r213=r212+6; //only col= 0 can read
    const int r221=r213+1; 
    const int r222=r221+6;
    const int r223=r222+6; //only col= 0 can read
//index to read old acceleration and force values, only col<4 can read this
    const int readOld=col+row*numBlocks*4+i1*4;

    __shared__ float Vals[6][4];
    __shared__ float temp1[6][6];
    __shared__ float temp2[6][6];
    __shared__ float A[6][6];

//indices to write solution
    const int wA11=col+row*numBlocks*8+i1*8;
    const int wF1c1=wA11+1;
    const int wA12=wF1c1+1;
    const int wF1c2=wA12+1;
    const int wA21=wF1c2+1;
    const int wF2c1=wA21+1;
    const int wA22=wF2c1+1;
    const int wF2c2=wA22+1;
    temp1[row][col]=0;
    temp2[row][col]=0;

    if(col<4)
    {
    Vals[row][col]=oldAs[readOld];
    }
    
    if(col==0)
    {
    temp1[row][col]=zs[r123];
    temp2[row][col]=zs[r213];
    }
    
    //temp1 is z23 body 1
    //temp2 is z13 body 2
    temp1[row][col]=temp1[row][col]-temp2[row][col];
    //temp1 is z123-z213
    syncthreads();
    if(col==0)
    {
    temp2[row][col]=Vals[row][3];
    }
    //temp2 is F2c2
    A[row][col]=zs[r212];
    syncthreads();
    Mat61Mult(A,temp2,temp2,row,col);
    temp2[row][col]=-1*temp2[row][col];
    //temp2 is -1*z212*F2c2
    temp2[row][col]=temp2[row][col]+temp1[row][col];
    //temp2 is -z212*F2c2+z123-z213
    syncthreads();
    temp1[row][col]=zs[r121];
    //temp1 is z121
    if(col==0)
    {
     A[row][col]=Vals[row][1];
    }
    //A is F1c1
    syncthreads();
    Mat61Mult(temp1,A,A,row,col);
    //A is z121*F1c1
    temp2[row][col]=A[row][col]+temp2[row][col];
    //temp2 is (z121*F1c1-z212*F2c2+z123-z213)
    if(row<5 && col<5)
    {
    A[row][col]=Xinv[rXinv];
    }
    make_W_pend(A,temp1,row,col);
    //A is now W
    Mat61Mult(A,temp2,A,row,col);
    A[row][col]=-1*A[row][col];
    //A is now -1*(W*((z121*A11)-(z212*F2c2+z123-z213)))
    if (col==0)
    {
    newAs[wF1c2]=A[row][col];//F1c2
    }
    
    temp1[row][col]=zs[r121];
    if(col==0)
    {
        temp2[row][col]=Vals[row][1];
    }
    syncthreads();
    Mat61Mult(temp1,temp2,temp1,row,col);
    temp2[row][col]=zs[r122];
    syncthreads();
    Mat61Mult(temp2,A,temp2,row,col);
    temp1[row][col]=temp1[row][col]+temp2[row][col];
    if(col==0)
    {
        temp2[row][col]=zs[r123];
    }
    syncthreads();
    temp2[row][col]=temp1[row][col]+temp2[row][col];//temp2 should hold A2 now
    
    if (col==0)
    {
        newAs[wA12]=temp2[row][col];//A12
    }

    
    A[row][col]=A[row][col]*-1;//Now A is the other F
    
    if (col==0)
    {
        newAs[wF2c1]=A[row][col];
    }
    
    temp1[row][col]=zs[r211];
    syncthreads();
    Mat61Mult(temp1,A,A,row,col);
    temp1[row][col]=zs[r212];
    if(col==0)
    {
    temp2[row][col]=Vals[row][3];
    }
    syncthreads();
    Mat61Mult(temp1,temp2,temp2,row,col);
    if (col==0)
    {
        temp1[row][col]=zs[r213];
    }
    temp1[row][col]=A[row][col]+temp2[row][col]+temp1[row][col];//temp1 now holds A1
    if(col==0)
    {
        newAs[wA21]=temp1[row][col];
        newAs[wA11]=Vals[row][0];
        newAs[wF1c1]=Vals[row][1];
        newAs[wA22]=Vals[row][2];
        newAs[wF2c2]=Vals[row][3];
    }
}

















__global__ void Assemble(float oldZs[],float newZs[],float Xinvs[], int nn)
{
    __shared__ float z1[6][6];
    __shared__ float z2[6][6];
    __shared__ float A[6][6];
    //Variable Declarations
    const int i1 = blockIdx.x;
    const int row = threadIdx.y;
    const int col = threadIdx.x;
    
    //indices to read the old zs
    const int r111=col+row*nn*52+i1*52;
    const int r112=r111+6;
    const int r113=r112+6; //only col = 0 can read
    const int r121=r113+1;
    const int r122=r121+6;
    const int r123=r122+6; //only col= 0 can read
    
    const int r211=r123+1;
    const int r212=r211+6;
    const int r213=r212+6; //only col= 0 can read
    const int r221=r213+1; 
    const int r222=r221+6;
    const int r223=r222+6; //only col= 0 can read
    
    int rowlen;
    rowlen = (nn)*26;
    //indices to write the new zs
    const int i11 = col+row*rowlen+i1*26;
    const int i12=i11+6;
    const int i13=i12+6; //only col=0 can write
    const int i21=i13+1;
    const int i22=i21+6;
    const int i23=i22+6; //only col = 0 can write
    
    //index to write Xinvs
    const int wXinv = col+row*(nn)*5+i1*5;//only write if col !=5 and row !=5
    
    /**************Get the intermediate quantity W***************/
    z1[row][col]=oldZs[r211];
    z2[row][col]=oldZs[r122];
    syncthreads();
    get_X_pend(z1,z2,A,row,col);
    
    invert_X(A,z1,row,col);
   
    if (col!=5 && row !=5)
    {
        Xinvs[wXinv]=A[row][col];
    }
    make_W_pend(A,z1,row,col);
    
   //A now holds W, so W can be used to find all newZs excpt z13 and z23
    z1[row][col]=oldZs[r112];//z12 for body 1
    z2[row][col]=oldZs[r212];//z12 for body 2
    syncthreads();
    Mat66Mult(A,z2,z2,row,col);
    Mat66Mult(z1,z2,z2,row,col);//z2 now holds the new z12
    newZs[i12]=z2[row][col];
    
    //next find z11
    z2[row][col]=oldZs[r121];//z21 for body 1
    syncthreads();
    Mat66Mult(A,z2,z2,row,col);
    Mat66Mult(z1,z2,z2,row,col);
    z1[row][col]=oldZs[r111];//z1 now holds z11 for body 1
    syncthreads();
    z2[row][col]=z1[row][col]-z2[row][col];//z2 now holds the new z11
    newZs[i11]=z2[row][col];
	
	//now z21
	z1[row][col]=oldZs[r221];
	z2[row][col]=oldZs[r121];
    syncthreads();
	Mat66Mult(A,z2,z2,row,col);
	Mat66Mult(z1,z2,z2,row,col);//z2 now holds the new z21
	newZs[i21]=z2[row][col];

    //now z22
	z2[row][col]=oldZs[r212];//z12 for body 2
    syncthreads();
	Mat66Mult(A,z2,z2,row,col);
	Mat66Mult(z1,z2,z2,row,col);
	z1[row][col]=oldZs[r222];
    syncthreads();
	z2[row][col]=z1[row][col]-z2[row][col];
	newZs[i22]=z2[row][col];

	//now find Y
	if(col==0)
	{
		z1[row][col]=oldZs[r123];
		z2[row][col]=oldZs[r213];
        
		
	}
    syncthreads();
    z1[row][col]=z1[row][col]-z2[row][col];
	Mat61Mult(A,z1,A,row,col);//col = 0 of A now holds Y

	//now find z13
	z2[row][col]=oldZs[r112];
    
	if (col==0)
	{
		z1[row][col]=oldZs[r113];
	}
    syncthreads();
	Mat61Mult(z2,A,z2,row,col);
	z1[row][col]=z1[row][col]-z2[row][col];
	if(col==0)
	{
		newZs[i13]=z1[row][col];
	}
	//now find z23
	z2[row][col]=oldZs[r221];
	if(col==0)
	{
		z1[row][col]=oldZs[r223];
	}
    syncthreads();
	Mat61Mult(z2,A,z2,row,col);
	z1[row][col]=z1[row][col]+z2[row][col];
	if(col==0)
	{
        newZs[i23]=z1[row][col];
	}
}

	
__device__ void get_W(float z1[6][6], float z2[6][6], float W[6][6],int row, int col)
{
    get_X_pend(z1,z2,W,row,col);
    invert_X(W,z1,row,col);
    make_W_pend(W,z1,row,col);
}

__device__ void get_X_pend(float z1[6][6], float z2 [6][6], float D[6][6],int row , int col)
{
    z1[row][col]=z1[row][col]+z2[row][col];
    
    set_D_pend(D,row,col);
    Mat66Mult(z1,D,z1,row,col);
    set_Dt_pend(D,row,col);
    Mat66Mult(D,z1,D,row,col);
}
__device__ void set_D_pend(float D[6][6],int row,int col)
{
    D[row][col]=0;
    if(row+col==0)
    {
        D[0][0]=1;
        D[1][1]=1;
        D[3][2]=1;
        D[4][3]=1;
        D[5][4]=1;
    }
}
__device__ void set_Dt_pend(float D[6][6],int row,int col)
{
    D[row][col]=0;
    if(row+col==0)
    {
        D[0][0]=1;
        D[1][1]=1;
        D[2][3]=1;
        D[3][4]=1;
        D[4][5]=1;
    }
}

__device__ void invert_X(float X[6][6], float A[6][6],int row, int col)
{
     int s;
     int j;
     float temp;
     float temp2;
     A[row][col]=0;
     if (row==col)
     {
          A[row][col]=1;
         for(s=0;s<5;s++)
         {
            temp=X[s][s];
            X[s][col]=(X[s][col]/temp);
            A[s][col]=(A[s][col]/temp);
            for(j=0;j<5;j++)
            {
                if(s != j)
                {
                     temp2=X[j][s];
                     X[j][col]=X[j][col]-(X[s][col]*temp2);
                     A[j][col]=A[j][col]-(A[s][col]*temp2);
                }
            }
          }
        }
        syncthreads();
        X[row][col]=A[row][col];
        syncthreads();
}
__device__ void make_W_pend(float Xinv[6][6], float D[6][6],int row , int col)                    
{
    if (row == 5 || col==5)
    {
         Xinv[row][col]=0;
    }
    
    set_Dt_pend(D,row,col);
    Mat66Mult(Xinv,D,Xinv,row,col);
    set_D_pend(D,row,col);
    Mat66Mult(D,Xinv,Xinv,row,col);
}















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

    int i = 0;
    while (i <= i1)
    {
        q+=state[i];
        w+=state[i+n];
        i++;
    }
   
    /*
    //or this, but only if less than 36 bodies:

    __shared__ float lcl_state[2*n];
    if (index < 2*n)
    {
        lcl_state[index]=state[index];
    }
    */
    //then add them up to get the right number

    ////////////////Inverse Mass Matrix and shifter setup//////////////////////////////////
    MSsetup(Minv,S,row,col,m[i1],I,Iindex,i1);
    
    
    
    ////////////////////DCM/////////////////////////////////////////
    //Dont need a dcm because it is planar, so for now to avoid memory 
    //usage, Im not using it.
    //DCM(q,C)
    
    ///////////////////z11/////////////////////////////////////////
    z[row][col]=0;
    r01(r,q,l[i1]);
 
    zaa(r,z,Minv,S,row,col);
    
    Zetas[i11]=z[row][col];
    //////////////////z12//////////////////////////////////////////
    //S is now S10, Minv is S10*Minv and z is z11
    z[row][col]=0;
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
    
    ////////////////z22///////////////////////////////////
    z[row][col]=0;
    MSsetup(Minv,S,row,col,m[i1],I,Iindex,i1);
    r02(r,q,l[i1]);
    zaa(r,z,Minv,S,row,col);
    Zetas[i22]=z[row][col];
    //////////////////z21/////////////////////////////////
    //S is now S20 Minv is S20*Minv and z is z22
    z[row][col]=0;
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
}
    
/**********************Zeta Calculation Functions**************/

//Function used to calculate z11 and z22
__device__ void zaa(float r[],float z[6][6],float Minv[6][6],float S[6][6], int row, int col)
{
     
     if ((row+col)==0)
     {
          makeS(S,r);
     }
   
   
     Mat66Mult(Minv,S,z,row,col);
     
     transpose(S,row,col);
     Mat66Mult(S,z,z,row,col);
     Mat66Mult(S,Minv,Minv,row,col);
}

//Function used to calculate z12 and z21
__device__ void zab(float r[],float z[6][6],float Minv[6][6],float S[6][6],int row , int col)
{    
    if (row+col==0)
     {
          makeS(S,r);
     }
    if (row>2 && col<3)
    {
         S[row][col]=0;
    }
    Mat66Mult(Minv,S,z,row,col);
}

/**********************Matrix Setup Functions******************/

//Creates a shifter matrix with the given position vector r
__device__ void makeS(float S[6][6],float r[])
{
     S[0][4]=-1*r[2];
     S[0][5]=r[1];
     S[1][3]=r[2];
     S[1][5]=-1*r[0];
     S[2][3]=-1*r[1];
     S[2][4]=r[0];
     
}


//Setup the inverse mass and shifter matrices.
__device__ void MSsetup(float Minv[6][6],float S[6][6],int row , int col,float m, float I[],int Iindex,int i1)
{
    Minv[row][col]=0;
    S[row][col]=0;
    if (row<3 && col==row)
    {
        Minv[row][col]=1/(I[Iindex]);
        Minv[row+3][col+3]=1/m;
    }
    if (col==row)
    {
        S[row][col]=1;
    }  
}

//set up the state dependant force vector Fa
__device__ void Fa(float S[6][6],float w, float r[], float m,int row, int col)
{
 float g = 9.81;
    if (col == 0)
    {
         S[row][col]=0;
    }
    if (row +col == 0)
    {
         S[3][0]=-m*w*w*r[0];
         S[4][0]=-g*m-m*w*w*r[1];
    }
}

//Function to find the r02 vector of the body
__device__ void r02(float r[], float angle,float L)
{
    r[0]=L*sin(angle)/2;
    r[1]=-1*L*cos(angle)/2;
    r[2]=0;
}

//Funtion to find the r01 position vector of the body
__device__ void r01(float r[], float angle,float L)
{
    r[0]=-1*L*sin(angle)/2;
    r[1]=L*cos(angle)/2;
    r[2]=0;
}

/**********************Matrix Operations*********************/

//Find the transpose of a 6x6 shifter matrix
__device__ void transpose(float S[6][6], int row, int col)
{
 if (row<3 && col>2)
 {
     S[col][row]=S[row][col];
     S[row][col]=0;
  }
}
//Multipy two 6x6 matrices, saving the result in C
__device__ void Mat66Mult(float A[6][6], float B[6][6], float C[6][6],int row, int col)
{
    float j = 0;
    int l=0;
    while (l<6)
    {
        j+=A[row][l]*B[l][col];
        l++;
    }
    syncthreads();
    C[row][col]=j;
    syncthreads();
}

//Multiply a 6x6 and a 6x1 matrix together.  The 6x1 is actually a 6x6 but 
//only the first collumn is used.
__device__ void Mat61Mult(float A[6][6], float B[6][6], float C[6][6],int row, int col)
{
    float j=0;
    int l =0;
    if (col==0)
    {
        while (l<6)
        {
            j+=A[row][l]*B[l][0];
            l++;
        }
        C[row][0]=j;
    }
    syncthreads();
}
/**********************Printing Functions**************************/

//Print a 6x6 matrix
__device__ void PrintMat66(float A[6][6], int nn )
{
syncthreads();
     int r;
     int c;
   if (blockIdx.x==nn)
   {
    if (threadIdx.y==0 && threadIdx.x==0)
    {    
         
         printf("\\n\\n");
         for( r=0;r<6;r++)
         {
              for(c=0;c<6;c++)
                {
                     printf("%f\\t",A[r][c]);
                }
                printf("\\n");
        }
        
      }
    }
    syncthreads();
}
"""
)
Initialize = mod.get_function("Initialize")
Assemble=mod.get_function("Assemble")
Disassemble=mod.get_function("Disassemble")


# In[97]:

# Initialization Function
#-------------------------
#This function prepares the system for the DCA algorithm 

def pycudaInitialize(b,x,n):
    
    x_gpu = x.astype(np.float32)
    m_gpu = np.zeros((n))
    l_gpu = np.zeros((n))
    I_gpu = np.zeros((3,3*n))
    for i in range (0,n):
        m_gpu[i]=b[i].m
        l_gpu[i]=b[i].l
        I_gpu[:,3*i:3*i+3]=b[i].I[:,:]
    m_gpu=m_gpu.astype(np.float32)
    l_gpu=l_gpu.astype(np.float32)
    I_gpu=I_gpu.astype(np.float32)
    gpu_out=np.zeros((6,26*n)).astype(np.float32)
    #time.sleep(.01)
    Initialize(cuda.In(x_gpu),cuda.In(m_gpu), cuda.In(l_gpu),cuda.In(I_gpu),cuda.Out(gpu_out),np.int32(n),block=(6,6,1),grid=(n,1,1))
    k=0
    index1 = 0
    for array in gpu_out:
        for element in array:
            if element < 1E-6:
                element = 0
    while k<n:
        b[k].z11=np.zeros((6,6))
        b[k].z12=np.zeros((6,6))
        b[k].z13=np.zeros((6))
        b[k].z21=np.zeros((6,6))
        b[k].z22=np.zeros((6,6))
        b[k].z23=np.zeros((6))
        b[k].z11[:,:]=gpu_out[:,index1:index1+6]
        index1 +=6
        b[k].z12[:,:]=gpu_out[:,index1:index1+6]
        index1 +=6
        b[k].z13[:]=gpu_out[:,index1]
        index1+=1
        b[k].z21[:,:]=gpu_out[:,index1:index1+6]
        index1+=6
        b[k].z22[:,:]=gpu_out[:,index1:index1+6]
        index1+=6
        b[k].z23[:]=gpu_out[:,index1]
        index1+=1
        k+=1


# In[98]:

def pycudaAssembly(bodies):
      
    odd = 0
    newbds=[]
    oldlen=np.int32(len(bodies))
    if (len(bodies)%2)==0:        
        newlen=math.trunc(oldlen/2)
        oldZs_gpu = np.zeros((6,26*oldlen))
        for k in range (0,newlen):
            newbds.append(MB.Body())
    else:
        odd=1
        oldZs_gpu = np.zeros((6,26*(oldlen-1)))
        newlen=math.trunc((oldlen/2)+1)
        for k in range (0,newlen):
            newbds.append(MB.Body())
    
    numBlocks=newlen-odd
    #Will have to mess with this for larger numbers of bodies and odd numbers
    newZs_gpu = np.zeros((6,26*(numBlocks))).astype(np.float32)
    Xinvs=np.zeros((5,5*numBlocks)).astype(np.float32)
    for i in range (0,numBlocks*2):
        k=26*i
        oldZs_gpu[:,k:k+6]=bodies[i].z11[:,:]
        oldZs_gpu[:,k+6:k+12]=bodies[i].z12[:,:]
        oldZs_gpu[:,k+12]=bodies[i].z13[:]
        oldZs_gpu[:,k+13:k+19]=bodies[i].z21[:,:]
        oldZs_gpu[:,k+19:k+25]=bodies[i].z22[:,:]
        oldZs_gpu[:,k+25]=bodies[i].z23[:]
    oldZs_gpu=oldZs_gpu.astype(np.float32)
    #time.sleep(.01)
    
    Assemble(cuda.In(oldZs_gpu),cuda.Out(newZs_gpu),cuda.Out(Xinvs),np.int32(numBlocks),block=(6,6,1),grid=(numBlocks,1,1))
    k=0
    index1 = 0
    index2= 0
    for array in newZs_gpu:
        for element in array:
            if element < 1E-6:
                element = 0
    for array in Xinvs:
        for element in array:
            if element < 1E-6:
                element = 0
                
    while k<numBlocks:
        newbds[k].z11=np.zeros((6,6))
        newbds[k].z12=np.zeros((6,6))
        newbds[k].z13=np.zeros((6))
        newbds[k].z21=np.zeros((6,6))
        newbds[k].z22=np.zeros((6,6))
        newbds[k].z23=np.zeros((6))
        newbds[k].Xinv=np.zeros((5,5))
        newbds[k].z11[:,:]=newZs_gpu[:,index1:index1+6]
        index1 +=6
        newbds[k].z12[:,:]=newZs_gpu[:,index1:index1+6]
        index1 +=6
        newbds[k].z13[:]=newZs_gpu[:,index1]
        index1+=1
        newbds[k].z21[:,:]=newZs_gpu[:,index1:index1+6]
        index1+=6
        newbds[k].z22[:,:]=newZs_gpu[:,index1:index1+6]
        index1+=6
        newbds[k].z23[:]=newZs_gpu[:,index1]
        index1+=1
        newbds[k].Xinv[:,:]=Xinvs[:,index2:index2+5]
        index2+=5
        k+=1
    if odd ==1:
        newbds[k].z11=np.zeros((6,6))
        newbds[k].z12=np.zeros((6,6))
        newbds[k].z13=np.zeros((6))
        newbds[k].z21=np.zeros((6,6))
        newbds[k].z22=np.zeros((6,6))
        newbds[k].z23=np.zeros((6))
        newbds[k].z11[:,:]=bodies[oldlen-1].z11[:,:]
        newbds[k].z12[:,:]=bodies[oldlen-1].z12[:,:]
        newbds[k].z13[:]=bodies[oldlen-1].z13[:]
        newbds[k].z21[:,:]=bodies[oldlen-1].z21[:,:]
        newbds[k].z22[:,:]=bodies[oldlen-1].z22[:,:]
        newbds[k].z23[:]=bodies[oldlen-1].z23[:]
        
    
    return newbds


# In[99]:

def pycudaDisassembly(oldAs,newbodies,oldbodies):
    
    odd = 0
    if (len(newbodies)%2)==0:  
       
        newlen=math.trunc(4*len(newbodies))
        newAs_gpu = np.zeros((6,newlen)).astype(np.float32)
    else:
        odd=1
       
        newlen=math.trunc(4*(len(newbodies)-1))
        newAs_gpu = np.zeros((6,newlen)).astype(np.float32)
        
        
    oldAs_gpu=np.zeros((6,len(oldAs)-4*odd))
    for k in range (0,len(oldAs)-4*odd):
        oldAs_gpu[:,k]=oldAs[k][:]
   
    numBlocks=np.int32(len(newbodies)-len(oldbodies))
    
   
    Zs_gpu = np.zeros((6,26*numBlocks*2))
    Xinvs=np.zeros((5,numBlocks*5))
    for i in range (0,numBlocks*2):
        k=26*i
        Zs_gpu[:,k:k+6]=newbodies[i].z11[:,:]
        Zs_gpu[:,k+6:k+12]=newbodies[i].z12[:,:]
        Zs_gpu[:,k+12]=newbodies[i].z13[:]
        Zs_gpu[:,k+13:k+19]=newbodies[i].z21[:,:]
        Zs_gpu[:,k+19:k+25]=newbodies[i].z22[:,:]
        Zs_gpu[:,k+25]=newbodies[i].z23[:]
        
    for k in range(0,numBlocks):
        i=5*k
        Xinvs[:,i:i+5]=oldbodies[k].Xinv[:,:]
    Xinvs=Xinvs.astype(np.float32)
    oldAs_gpu=oldAs_gpu.astype(np.float32)
    Zs_gpu=Zs_gpu.astype(np.float32)
    #time.sleep(.01)
    Disassemble(cuda.In(Xinvs),cuda.In(Zs_gpu),cuda.In(oldAs_gpu),cuda.Out(newAs_gpu),numBlocks,block=(6,6,1),grid=(int(numBlocks),1,1))
    
    for array in newAs_gpu:
        for element in array:
            if element < 1E-6:
                element = 0
    newsol=[]           
    for k in range (0,newlen+4*odd):
        newsol.append(np.zeros((6)))
    k=0
    while k<newlen:
        newsol[k][:]=newAs_gpu[:,k]
        k+=1
    if odd ==1:
        newsol[k][:]=oldAs[-4][:]
        k+=1
        newsol[k][:]=oldAs[-3][:]
        k+=1
        newsol[k][:]=oldAs[-2][:]
        k+=1
        newsol[k][:]=oldAs[-1][:]
    
    return newsol


# In[100]:

# Recursive DCA Algorithm
#------------------------
#This function applies the DCA Algorithm to a compound pendulum 
#system, regardless of the number of links.  Recursion is 
#used to assemble and disassemble the system.

def recursiveDCA(n,i,bodies,joints,BC1,BC2,xinits): 
    
    #This if statements is the recursive case of this function.
    #Inside this if statement the current list of bodies is
    #assembled into a smaller list and passed back to this function.
    if len(bodies) !=1: #if there are still more than one bodies in the system
        odd = 0
        
        #check if there are an odd number of bodies and
        #create a new list called "newbds" containing the correct number of
        #new bodies
        if (len(bodies)%2)!=0:
            odd=1
           
            
        newbds=pycudaAssembly(bodies)
        
        #Find the new joints corresponding to the new bodies
        newjs=[]
        for k in range (0,len(newbds)):
            newjs.append(MB.Joint())
        
        
        #loop that brings along the P and D matrices from the previous joints
        for j in range(0,len(newjs)):
            if j==len(newjs)-1 and odd ==1:
                newjs[j].D=joints[len(joints)-1].D
                newjs[j].P=joints[len(joints)-1].P
               
            else:
                newjs[j].D=joints[2*j].D
                newjs[j].P=joints[2*j].P
               
                
                
        #This is the recursive call.  This will return a list of the form:
        #[A11,F1c1,A12,F1c2,A21,F2c1,...,An2,Fnc2] where the Axy's and Fxcy's 
        #correspond to the accellerations and constraint forces in the 
        #xth body on the yth handle. The function is given the original 
        #number of bodies, the list of new bodies, the list of new joints,
        #the boundary conditions,and the state variables at that timestep.
        #At this point this function will repeat itself in a loop until there is
        #one body left

        sol=recursiveDCA(n,i+1,newbds,newjs,BC1,BC2,xinits)
        
        newsol=pycudaDisassembly(sol,bodies,newbds)
        #If this is the 0th level of recursion, delete all the forces 
        #so only the accellerations are returned
        if i ==0:
            for j in range (1,2*n+1):
                del newsol[j]
                  
            return newsol
        
        #If this is not the 0th level of recursion, return the solution
        #with all of the forces included
        else:
            return newsol
    
    #This is the base case of this recursive function.  This is where
    #the recursion ends and reassembly begins
    elif BC1==2 and BC2==1:
        
        #Fc2 is zero because that end is free
        Fc2=np.zeros((6))
        
        #Because the first joint is a pin, it will
        #have no translational motion, and support no 
        #moments, so mat is used to solve for the two translational
        #forces of Fc1
        mat=np.zeros((3,4))
        mat[:,:3]=bodies[0].z11[3:,3:]
        mat[:,3]=-1*bodies[0].z13[3:]
        
        #Loop to put a matrix in reduced row echelon form
        for s in range(0,3):
            mat[s,:]=mat[s,:]/mat[s,s]
            for j in range(0,3):
                if s !=j:
                    mat[j,:]=mat[j,:]-(mat[s,:]*mat[j,s])
    
        Fc1=np.zeros_like(Fc2)
        Fc1[3:]=mat[:,3]

        #solve for the A's given Fc2=0
        A1=np.dot(bodies[0].z11,Fc1)+bodies[0].z13
        A2=np.dot(bodies[0].z21,Fc1)+bodies[0].z23
        sol=[A1,Fc1,A2,Fc2]
        
        #return the forces and accellerations at the end joints
        #and begin disassembly
        return sol  


# In[101]:

# Helper Function
#-------------------------
#This function provides a buffer between the dca algorithm and the 
#integrator call. 
def funct(state,t,n,i,bbs,jjs,BC1,BC2):

    #Initialize the system
    pycudaInitialize(bbs,state,n)
    #Call the Recursive DCA Algorithm
    #This returns a list of the form:
    #[A11,A12,A21,A22,...,An1,An2]
    #where Axy corresponds to the acceleration
    #of the yth handle of the xth body
    ll=recursiveDCA(n,0,bbs,jjs,BC1,BC2,state)
    
    #loop to fill d_dt with the acceleration values
    d_dt=np.zeros((2*n))
    for j in range(0,n):
        if j == 0:
            A1=ll.pop(0)
            d_dt[j+n]=np.dot(np.transpose(jjs[j].P),A1)
        else:
            A2= ll.pop(0)
            A1=ll.pop(0)
            d_dt[j+n]=np.dot(np.transpose(jjs[j].P),(A1-A2))
    
    #add the velocities to d_dt and return to the integrator
    d_dt[:n]=state[n:]
   
    return d_dt 


# In[102]:

# Integration
#-------------

#odeint is the numerical integrator used
yy=odeint(funct,initialconditions,Time,(n,0,bs,js,2,1))


# In[103]:

# Energy Calculation
#--------------------
#The energy of the system is calculated and plotted

energy=MBF.PendEnergy(yy,bs)
KE=energy[:,0]
PE=energy[:,1]
TE=energy[:,2]

plt.plot(Time,TE-TE[0])
plt.xlabel("Time [s]")
plt.ylabel("energy")
plt.title("System Energy")
plt.show()

plt.plot(Time,PE,Time,KE)
plt.xlabel("Time[s]")
plt.ylabel("energy")
plt.title("Kinetic and Potential Energy")
plt.show()


# In[104]:

# Solution Plot
#--------------------
plt.plot(Time,yy[:,:n])
plt.xlabel("Time [s]")
plt.ylabel("Generalized Coordinates [Rad]")
plt.title("System Response")

plt.show()

plt.plot(Time,yy[:,n:])
plt.xlabel(("Time[s]"))
plt.ylabel(("Generalized Speeds [Rad/s]"))
plt.title("System Response")

plt.show()


# In[105]:

# Recursive DCA Algorithm
#------------------------
#This function applies the DCA Algorithm to a compound pendulum 
#system, regardless of the number of links.  Recursion is 
#used to assemble and disassemble the system.

def recursiveDCA2(n,i,bodies,joints,BC1,BC2,xinits): 
    
    #This if statements is the recursive case of this function.
    #Inside this if statement the current list of bodies is
    #assembled into a smaller list and passed back to this function.
    if len(bodies) !=1: #if there are still more than one bodies in the system
        j=0
        odd = 0
        
        #check if there are an odd number of bodies and
        #create a new list called "newbds" containing the correct number of
        #new bodies
        newbds=[]
        if (len(bodies)%2)==0:
            for k in range (0,math.trunc(len(bodies)/2)):
                newbds.append(MB.Body())
        else:
            odd=1
            for k in range (0,math.trunc(len(bodies)/2)+1):
                newbds.append(MB.Body())
            
        #Loop through all of the new bodies, assembling the correct two
        #"old" bodies to create the new.
        while j < len(newbds):
            
            #if there are an odd # of bodies and we're at the last 
            #body in newbds, add the last body to newbds
            if j == len(newbds)-1 and odd == 1:
                newbds[j].z11=np.zeros((6,6))
                newbds[j].z12=np.zeros((6,6))
                newbds[j].z13=np.zeros((6))
                newbds[j].z21=np.zeros((6,6))
                newbds[j].z22=np.zeros((6,6))
                newbds[j].z23=np.zeros((6))
                newbds[j].z11[:,:]=bodies[2*j].z11[:,:]
                newbds[j].z12[:,:]=bodies[2*j].z12[:,:]
                newbds[j].z21[:,:]=bodies[2*j].z21[:,:]
                newbds[j].z22[:,:]=bodies[2*j].z22[:,:]
                newbds[j].z13[:]=bodies[2*j].z13[:]
                newbds[j].z23[:]=bodies[2*j].z23[:]
                newbds[j].m=bodies[2*j].m
                
            #Otherwise, calculate the new zetas for the newly assembled bodies
            #according to the formulae
            else:
                newbds[j].m=bodies[2*j].m+bodies[2*j+1].m         
                #X, Q, and Y  are  intermediate quantities
                newbds[j].X=np.dot(np.transpose(joints[2*j+1].D),np.dot(bodies[2*j+1].z11+bodies[2*j].z22,joints[2*j+1].D))
                newbds[j].Xinv=np.linalg.inv(newbds[j].X)
                newbds[j].W=np.dot(joints[2*j+1].D,np.dot(newbds[j].Xinv,np.transpose(joints[2*j+1].D)))
                newbds[j].Y=np.dot(newbds[j].W,bodies[2*j].z23-bodies[2*j+1].z13)#ommitted pdot*u
                #assemble the bodies based on the formulas
                newbds[j].z11=bodies[2*j].z11-np.dot(bodies[2*j].z12,np.dot(newbds[j].W,bodies[2*j].z21))
                newbds[j].z12=np.dot(bodies[2*j].z12,np.dot(newbds[j].W,bodies[2*j+1].z12))
                newbds[j].z21=np.dot(bodies[2*j+1].z21,np.dot(newbds[j].W,bodies[2*j].z21))
                newbds[j].z22=bodies[2*j+1].z22-np.dot(bodies[2*j+1].z21,np.dot(newbds[j].W,bodies[2*j+1].z12))
                newbds[j].z13=bodies[2*j].z13-np.dot(bodies[2*j].z12,newbds[j].Y)
                newbds[j].z23=bodies[2*j+1].z23+np.dot(bodies[2*j+1].z21,newbds[j].Y)
            j=j+1
            
        
        #Find the new joints corresponding to the new bodies
        newjs=[]
        for k in range (0,len(newbds)):
            newjs.append(MB.Joint())
        
        
        #loop that brings along the P and D matrices from the previous joints
        for j in range(0,len(newjs)):
            if j==len(newjs)-1 and odd ==1:
                newjs[j].D=joints[len(joints)-1].D
                newjs[j].P=joints[len(joints)-1].P
               
            else:
                newjs[j].D=joints[2*j].D
                newjs[j].P=joints[2*j].P
               
                
                
        #This is the recursive call.  This will return a list of the form:
        #[A11,F1c1,A12,F1c2,A21,F2c1,...,An2,Fnc2] where the Axy's and Fxcy's 
        #correspond to the accellerations and constraint forces in the 
        #xth body on the yth handle. The function is given the original 
        #number of bodies, the list of new bodies, the list of new joints,
        #the boundary conditions,and the state variables at that timestep.
        #At this point this function will repeat itself in a loop until there is
        #one body left

        sol=recursiveDCA(n,i+1,newbds,newjs,BC1,BC2,xinits)
       
        
        #Forces and Accelerations at the new joints are found,
        #so these values can now be found for the old joints.
        
        #newsol will contain the new solution
        newsol=[];
        newsol.append(sol[0])
        newsol.append(sol[1])
        
       
        
        #In this loop I start with the first body and find the force and 
        #acceleration on its other handle. These values are added to the new
        #solution list, along with the next values in sol.
        flag=0
        for j in range(0,len(newbds)):
            
            #Don't enter this if there are an odd number of bodies and it 
            #is the last time through the loop, otherwise enter.
            if not(odd==1 and j==len(newbds)-1):
                
                #index counter
                k=len(newsol)-2
                
                #The force in the joint between two assembled bodies
                #F=np.dot(np.linalg.inv(bodies[2*j].z12),newsol[k]-bodies[2*j].z13-np.dot(bodies[2*j].z11,newsol[k+1]))
                F=-1*(np.dot(joints[j+1].D,np.dot(newbds[j].Xinv,np.dot(np.transpose(joints[j+1].D),(np.dot(bodies[2*j].z21,newsol[k+1])-np.dot(bodies[2*j+1].z12,sol[4*j+3])+bodies[2*j].z23-bodies[2*j+1].z13)))))
                
                #A2 is the acceleration of handle 2 of body k
                A2=np.dot(bodies[2*j].z21,newsol[k+1])+np.dot(bodies[2*j].z22,F)+bodies[2*j].z23
               
                #A1 is the acceleration of handle 1 of body k+1
                A1=np.dot(bodies[2*j+1].z11,-1*F)+np.dot(bodies[2*j+1].z12,sol[4*j+3])+bodies[2*j+1].z13
                 
                #Add the newly found values to new solution list,
                #maintiaining the format described above.
                newsol.append(A2)
                newsol.append(F)
                newsol.append(A1)
                newsol.append(-1*F)
                
                #Add the next two values in sol to the new solution list
                newsol.append(sol[4*j+2])
                newsol.append(sol[4*j+3])
            
            #if there are an odd number of bodies and this is the last
            #time through the loop, append the last force and accelleration 
            #in the old solution to the new one
            else:
                newsol.append(sol[4*j+2])
                newsol.append(sol[4*j+3])
                flag=1
            
            #If this is not the last time through the loop, append the 
            #next two values in the old solution to the new solution
            if j!=len(newbds)-1  and flag !=1:   
                newsol.append(sol[4*j+4])
                newsol.append(sol[4*j+5])
        
        #If this is the 0th level of recursion, delete all the forces 
        #so only the accellerations are returned
        if i ==0:
            for j in range (1,2*n+1):
                del newsol[j]
                  
            return newsol
        
        #If this is not the 0th level of recursion, return the solution
        #with all of the forces included
        else:
            return newsol
    
    #This is the base case of this recursive function.  This is where
    #the recursion ends and reassembly begins
    elif BC1==2 and BC2==1:
        
        #Fc2 is zero because that end is free
        Fc2=np.zeros((6))
        
        #Because the first joint is a pin, it will
        #have no translational motion, and support no 
        #moments, so mat is used to solve for the two translational
        #forces of Fc1
        mat=np.zeros((3,4))
        mat[:,:3]=bodies[0].z11[3:,3:]
        mat[:,3]=-1*bodies[0].z13[3:]
        
        #Loop to put a matrix in reduced row echelon form
        for s in range(0,3):
            mat[s,:]=mat[s,:]/mat[s,s]
            for j in range(0,3):
                if s !=j:
                    mat[j,:]=mat[j,:]-(mat[s,:]*mat[j,s])
    
        Fc1=np.zeros_like(Fc2)
        Fc1[3:]=mat[:,3]

        #solve for the A's given Fc2=0
        A1=np.dot(bodies[0].z11,Fc1)+bodies[0].z13
        A2=np.dot(bodies[0].z21,Fc1)+bodies[0].z23
        sol=[A1,Fc1,A2,Fc2]
        
        #return the forces and accellerations at the end joints
        #and begin disassembly
        return sol  


# In[106]:

# Initialization Function
#-------------------------
#This function prepares the system for the DCA algorithm 

def initialize(b,x,n):
    val1=0
    val2=np.zeros((3))
    val3=np.zeros((3))
    for k in range(0,n):
        val1=val1+x[k]
        b[k].C0=MBF.simdcm(val1,np.array((0,0,1)))
        val2=val2+np.array((0,0,x[k+n]))
        b[k].w=val2
        b[k].w1=np.array((0,0,x[k+n]))
        b[k].rs(b[k].l/2.0,np.array((0,-1,0)))
        b[k].Vcm=val3+np.cross(b[k].w,b[k].r10)
        val3=val3+np.cross(b[k].w,b[k].r12)
        b[k].Inertia2()#transformed inertia matrix
        b[k].Mmatrix()#Mass matrix and its inverse
        b[k].Forces(np.array((0,-1,0)))
        b[k].shifters()#S matrices
        b[k].zs()#zeta values


# In[107]:

# Helper Function
#-------------------------
#This function provides a buffer between the dca algorithm and the 
#integrator call. 
def funct2(state,t,n,i,bbs,jjs,BC1,BC2):

    #Initialize the system
    initialize(bbs,state,n)
    
    #Call the Recursive DCA Algorithm
    #This returns a list of the form:
    #[A11,A12,A21,A22,...,An1,An2]
    #where Axy corresponds to the acceleration
    #of the yth handle of the xth body
    ll=recursiveDCA2(n,0,bbs,jjs,BC1,BC2,state)
    
    #loop to fill d_dt with the acceleration values
    d_dt=np.zeros((2*n))
    for j in range(0,n):
        if j == 0:
            A1=ll.pop(0)
            d_dt[j+n]=np.dot(np.transpose(jjs[j].P),A1)
        else:
            A2= ll.pop(0)
            A1=ll.pop(0)
            d_dt[j+n]=np.dot(np.transpose(jjs[j].P),(A1-A2))
    
    #add the velocities to d_dt and return to the integrator
    d_dt[:n]=state[n:]
    return d_dt 


# In[108]:

# Integration
#-------------

#odeint is the numerical integrator used
yy2=odeint(funct2,initialconditions,Time,(n,0,bs,js,2,1))
print(yy2.astype(np.float16)-yy.astype(np.float16))


# In[109]:

# Energy Calculation
#--------------------
#The energy of the system is calculated and plotted

energy=MBF.PendEnergy(yy2,bs)
KE=energy[:,0]
PE=energy[:,1]
TE=energy[:,2]

plt.plot(Time,TE-TE[0])
plt.xlabel("Time [s]")
plt.ylabel("energy")
plt.title("System Energy")
plt.show()

plt.plot(Time,PE,Time,KE)
plt.xlabel("Time[s]")
plt.ylabel("energy")
plt.title("Kinetic and Potential Energy")
plt.show()


# In[110]:

# Solution Plot
#--------------------
plt.plot(Time,yy2[:,:n])
plt.xlabel("Time [s]")
plt.ylabel("Generalized Coordinates [Rad]")
plt.title("System Response")

plt.show()

plt.plot(Time,yy2[:,n:])
plt.xlabel(("Time[s]"))
plt.ylabel(("Generalized Speeds [Rad/s]"))
plt.title("System Response")

plt.show()

