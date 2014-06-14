__device__ void set_Dt_pend(float D[6][6],int row,int col);
__device__ void make_W_pend(float Xinv[6][6], float D[6][6],int row , int col);                    
__device__ void get_X_pend(float z1[6][6], float z2 [6][6], float D[6][6],int row , int col);
__device__ void set_D_pend(float D[6][6],int row,int col);
__device__ void get_W(float z1[6][6], float z2[6][6], float W[6][6],int row, int col);
__device__ void invert_X(float X[6][6], float A[6][6],int row, int col);
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

#include "../d_code/deviceDisassemble.h"
#include "../d_code/deviceAssemble.h"
#include "../d_code/deviceInitialize.h"
#include "../d_code/deviceFuncts.h"
