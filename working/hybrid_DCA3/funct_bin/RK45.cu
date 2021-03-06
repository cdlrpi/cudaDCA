//	RK45.cu
//
//This file contains the function that performs Runge Kutta 45 integration using the DCA algorithm

//Included Files
#include "classes.h"
#include <iostream>

//Function Prototypes
//	Functions found in Functs.cu
void arycpy(double A[],double B[],int n);
void arycpy2(double A[],double B[],int n);
void arycpy3(double A[],double B[],int n);
void set_up(double A[], double B[], double C[], int n , double h);
void get_final_q(double x[], double h, double a[], double b[], double c[], double d[], double R[], int n);
void get_final_qdot(double x[], double h, double a[], double b[], double c[], double d[], double R[], int n);
void DCAhelp(double state[], double m[], double l[], double I[],int n, double Y[],int cut_off, double Zs[], float times[],int reps);


//RK_45
//	Function used to perform Runge-Kutta 45 integration using the DCA algorithm.
//	This numerical integrator uses a fixed time step.
//		state is an array of the conditions at that timestep
//		step is the length of one timestep
//		n is the number of bodies
//		bs is a list of bodies
//		js is a list of joints
//		Y is the array where the conditions at the next timestep will be stored
void RK_45(double state[], double step, int n,double m[], double l[], double II[],double Y[],int cut_off, double Zs[],float times[], int reps)
{
	//Comments have not yet been completed in this file because I do not know how it works
	//Variable Declarations
	double *q = new double[n];
	double *qdot = new double[n];
	double *q1 = new double[n];
	double *qdot1 = new double[n];
	double *q2 = new double[n];
	double *qdot2 = new double[n];
	double *q3 = new double[n];
	double *qdot3 = new double[n];
	double *q4 = new double[n];
	double *qdot4 = new double[n];
	double *qddot1 = new double[n];
	double *qddot2 = new double[n];
	double *qddot3 = new double[n];
	double *qddot4 = new double[n];
	double *tState = new double[2*n];
	double *Y1 = new double[2*n];
	double *Y2= new double[2*n];
	double *Y3=new double[2*n];
	double *Y4=new double[2*n];
	int i = 0;
	double h = step/2.0;
	
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
	DCAhelp(tState,m,l,II,n,Y1,cut_off,Zs,times,reps);

	arycpy3(qddot1,Y1,n);
	
/*
		
	set_up(q,qdot1,q2,n,h);
	set_up(qdot,qddot1,qdot2,n,h);
	arycpy(tState,q2,n);
	arycpy2(tState,qdot2,n);
	DCAhelp(tState,m,l,II,n,Y2,cut_off,Zs);
	arycpy3(qddot2,Y2,n);
	set_up(q,qdot2,q3,n,h);
	set_up(qdot,qddot2,qdot3,n,h);
	arycpy(tState,q3,n);
	arycpy2(tState,qdot3,n);
	DCAhelp(tState,m,l,II,n,Y3,cut_off,Zs);
	arycpy3(qddot3,Y3,n);
	set_up(q,qdot3,q4,n,h*2);
	set_up(qdot,qddot3,qdot4,n,h*2);
	arycpy(tState,q4,n);
	arycpy2(tState,qdot4,n);
	DCAhelp(tState,m,l,II,n,Y4,cut_off,Zs);
	arycpy3(qddot4,Y4,n);
*/

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

void RK_452(double state[], double step, int n,double m[], double l[], double II[],double Y[],int cut_off, double Zs[])
{
	//Comments have not yet been completed in this file because I do not know how it works
	//Variable Declarations
	double *q = new double[n];
	double *qdot = new double[n];
	double *q1 = new double[n];
	double *qdot1 = new double[n];
	double *q2 = new double[n];
	double *qdot2 = new double[n];
	double *q3 = new double[n];
	double *qdot3 = new double[n];
	double *q4 = new double[n];
	double *qdot4 = new double[n];
	double *qddot1 = new double[n];
	double *qddot2 = new double[n];
	double *qddot3 = new double[n];
	double *qddot4 = new double[n];
	double *tState = new double[2*n];
	double *Y1 = new double[2*n];
	double *Y2= new double[2*n];
	double *Y3=new double[2*n];
	double *Y4=new double[2*n];
	int i = 0;
	double h = step/2.0;
	
	while ( i < n)
	{
		q[i]=state[i];
		qdot[i]=state[i+n];
		i++;
	}
	
	float times[3];
	int reps;
	arycpy(q1,q,n);
	arycpy(qdot1,qdot,n);
	arycpy(tState,q1,n);
	arycpy2(tState,qdot1,n);
	DCAhelp(tState,m,l,II,n,Y1,cut_off,Zs,times,reps);

	arycpy3(qddot1,Y1,n);
	

		
	set_up(q,qdot1,q2,n,h);
	set_up(qdot,qddot1,qdot2,n,h);
	arycpy(tState,q2,n);
	arycpy2(tState,qdot2,n);
	DCAhelp(tState,m,l,II,n,Y2,cut_off,Zs,times,reps);
	arycpy3(qddot2,Y2,n);
	set_up(q,qdot2,q3,n,h);
	set_up(qdot,qddot2,qdot3,n,h);
	arycpy(tState,q3,n);
	arycpy2(tState,qdot3,n);
	DCAhelp(tState,m,l,II,n,Y3,cut_off,Zs,times,reps);
	arycpy3(qddot3,Y3,n);
	set_up(q,qdot3,q4,n,h*2);
	set_up(qdot,qddot3,qdot4,n,h*2);
	arycpy(tState,q4,n);
	arycpy2(tState,qdot4,n);
	DCAhelp(tState,m,l,II,n,Y4,cut_off,Zs,times,reps);
	arycpy3(qddot4,Y4,n);


	get_final_q(q,h,qdot1,qdot2,qdot3,qdot4,Y, n);
	get_final_qdot(qdot,h,qddot1,qddot2,qddot3,qddot4,Y, n);
	
	std::cout<<Y[2]<<'\t'<<Y[3]<<std::endl;
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
