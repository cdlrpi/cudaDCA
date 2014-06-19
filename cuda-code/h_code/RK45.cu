
#include "classes.h"
#include <iostream>
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
void DCAhelp(float state[], InitBody *bs, Joint *js,int n, float Y[]);
void printa(float A[], int n);
//Function to perform Runge-Kutta 45 integration
void RK_45(float state[], float step, int n, InitBody *bs, Joint *js,float Y[])
{
	float *q = new float[n];
	float *qdot = new float[n];
	float *q1 = new float[n];
	float *qdot1 = new float[n];
	float *q2 = new float[n];
	float *qdot2 = new float[n];
	float *q3 = new float[n];
	float *qdot3 = new float[n];
	float *q4 = new float[n];
	float *qdot4 = new float[n];
	float *qddot1 = new float[n];
	float *qddot2 = new float[n];
	float *qddot3 = new float[n];
	float *qddot4 = new float[n];
	float *tState = new float[2*n];
	float *Y1 = (float*) malloc(sizeof(float)*2*n);
	float *Y2= new float[2*n];
	float *Y3=new float[2*n];
	float *Y4=new float[2*n];
	int i = 0;
	float h = step/2.0;
	std::cout<<h<<std::endl;
	while ( i < n)
	{
		q[i]=state[i];
		qdot[i]=state[i+n];
		i++;
	}
	
	arycpy(q1,q,n);
	arycpy(qdot1,qdot,n);
	arycpy(tState,q1,n);
	arycpy2(tState,qdot1,n);
	DCAhelp(tState,bs,js,n,Y1);
	arycpy3(qddot1,Y1,n);
	//printa(Y1,2*n);
	set_up(q,qdot1,q2,n,h);
	set_up(qdot,qddot1,qdot2,n,h);
	arycpy(tState,q2,n);
	arycpy2(tState,qdot2,n);
	std::cout<<"this"<<std::endl;
	//printa(tState,2*n);
	DCAhelp(tState,bs,js,n,Y2);
	//printa(Y2,2*n);	
	arycpy3(qddot2,Y2,n);
	

	
	set_up(q,qdot2,q3,n,h);
	set_up(qdot,qddot2,qdot3,n,h);
	arycpy(tState,q3,n);
	arycpy2(tState,qdot3,n);
	//printa(tState,2*n);	
	DCAhelp(tState,bs,js,n,Y3);
	arycpy3(qddot3,Y3,n);
	//printa(tState,2*n);


	set_up(q,qdot3,q4,n,h*2);
	set_up(qdot,qddot3,qdot4,n,h*2);
	arycpy(tState,q4,n);
	arycpy2(tState,qdot4,n);
	//printa(tState, 2*n);
	DCAhelp(tState,bs,js,n,Y4);
	arycpy3(qddot4,Y4,n);


	get_final_q(q,h,qdot1,qdot2,qdot3,qdot4,Y, n);
	get_final_qdot(qdot,h,qddot1,qddot2,qddot3,qddot4,Y, n);
	
	delete[] q;
	delete[] qdot;
	delete[] q1;
	delete[] qdot1;
	delete[] q2;
	delete[] qdot2;
	delete[] q3;
	delete[] qdot3;
	delete[] q4;
	delete[] qdot4;
	delete[] qddot1;
	delete[] qddot2;
	delete[] qddot3;
	delete[] qddot4;
	delete[] tState;
	delete[] Y1;
	delete[] Y2;
	delete[] Y3;
	delete[] Y4;
}

