__device__ void set_Dt_pend(double D[6][6],int row,int col);
__device__ void make_W_pend(double Xinv[6][6], double D[6][6],int row , int col);                    
__device__ void get_X_pend(double z1[6][6], double z2 [6][6], double D[6][6],int row , int col);
__device__ void set_D_pend(double D[6][6],int row,int col);
__device__ void get_W(double z1[6][6], double z2[6][6], double W[6][6],int row, int col);
__device__ void invert_X(double X[6][6], double A[6][6],int row, int col);
__device__ void MSsetup(double Minv[6][6],double S[6][6],int row, int col, double m, double I[],int Iindex,int i1);
__device__ void Fa(double S[6][6],double w, double r[], double m,int row, int col);
__device__ void transpose(double S[6][6],int row, int col);
__device__ void makeS(double S[6][6],double *r);
__device__ void zab(double r[],double z[6][6],double Minv[6][6],double S[6][6], int row, int col);
__device__ void zaa(double r[],double z[6][6],double Minv[6][6],double S[6][6], int row, int col);
__device__ void r02(double r[], double angle,double L);
__device__ void r01(double r[], double angle,double L);
__device__ void DCM(double angle, double **B);
__device__ void Mat61Mult(double A[6][6], double B[6][6], double C[6][6],int row, int col);
__device__ void Mat66Mult(double A[6][6], double B[6][6], double C[6][6],int row, int col);
__device__ void printit(double A[6][6]);

#include "../d_code/deviceDisassemble.h"
#include "../d_code/deviceAssemble.h"
#include "../d_code/deviceInitialize.h"
#include "../d_code/deviceFuncts.h"
