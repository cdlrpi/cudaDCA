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
	__syncthreads();
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
	__syncthreads();
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
        __syncthreads();
        X[row][col]=A[row][col];
        __syncthreads();
}
__device__ void make_W_pend(float Xinv[6][6], float D[6][6],int row , int col)                    
{
    if (row == 5 || col==5)
    {
         Xinv[row][col]=0;
    }
    	__syncthreads();
    set_Dt_pend(D,row,col);
    Mat66Mult(Xinv,D,Xinv,row,col);
    set_D_pend(D,row,col);
    Mat66Mult(D,Xinv,Xinv,row,col);
}
/**********************Zeta Calculation Functions**************/

//Function used to calculate z11 and z22
__device__ void zaa(float r[],float z[6][6],float Minv[6][6],float S[6][6], int row, int col)
{
     
     if ((row+col)==0)
     {
          makeS(S,r);
     }
   	__syncthreads();
   
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
	__syncthreads();
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
		__syncthreads();  
}

//set up the state dependant force vector Fa
__device__ void Fa(float S[6][6],float w, float r[], float m,int row, int col)
{
 float g = 9.81;
    if (col == 0)
    {
         S[row][col]=0;
    }
	__syncthreads();
    if (row +col == 0)
    {
         S[3][0]=-m*w*w*r[0];
         S[4][0]=-g*m-m*w*w*r[1];
    }
	__syncthreads();
}

//Function to find the r02 vector of the body
__device__ void r02(float r[], float angle,float L)
{
    r[0]=L*sin(angle)/2;
    r[1]=-1*L*cos(angle)/2;
    r[2]=0;
	__syncthreads();
}

//Funtion to find the r01 position vector of the body
__device__ void r01(float r[], float angle,float L)
{
    r[0]=-1*L*sin(angle)/2;
    r[1]=L*cos(angle)/2;
    r[2]=0;
	__syncthreads();
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
	__syncthreads();
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
    __syncthreads();
    C[row][col]=j;
    __syncthreads();
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
    __syncthreads();
}

