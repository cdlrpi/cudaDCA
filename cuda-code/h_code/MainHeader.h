#include "deviceHeader.h"
#include <iostream>
#include <math.h>
#include "classes.h"
#include "HelperFuncts.h"

//Function Prototypes
void pend_init(InitBody *bs,int n,float mass, float length);
void full_drop(float x[],int n);
void Time_init(float time[], float step, float final);
void RK_45(float state[], float step, int n, InitBody *bs, Joint *js,float Y[]);
void arycpy(float A[],float B[],int n);
void arycpy2(float A[],float B[],int n);
void arycpy3(float A[],float B[],int n);
void set_up(float A[], float B[], float C[], int n , float h);
void get_final_q(float x[], float h, float a[], float b[], float c[], float d[], float R[], int n);
void get_final_qdot(float x[], float h, float a[], float b[], float c[], float d[], float R[], int n);
__global__ void Initialize(float state[],float m[], float l[],float I[],float Zetas[],int n);
void DCAhelp(float state[], InitBody *bs, Joint *js,int n, float Y[]);
void CudaInitialize(InitBody* oldbs, float x[], int n, Body *newbs);
void printm(float A[6][6]);
void printa(float A[], int x);
__global__ void Assemble(float oldZs[],float newZs[],float Xinvs[], int nn);
void cudaAssemble(Body *bodies, int num, Body *newbds, int odd, int newlen);
void RecDCA(Body *bodies, int n, int i, float AF[]);
__global__ void Disassemble(float Xinv[],float Zs[],float oldAF[], float newAF[], int numBlocks);
void cudaDisassemble(float OldAF[], Body *morebds, Body *lessbds, int odd, int morelen, int lesslen, float AF[]);
void solve_BCs(Body *bodies, float AF[]);
void matmult61(float A[6][6], float B[6], float C[6]);
