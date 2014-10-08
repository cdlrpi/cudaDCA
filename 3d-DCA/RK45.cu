//	RK45.cu
//
//This file contains the function that performs Runge Kutta 45 integration using the DCA algorithm

//Included Files

#include <iostream>

//Function Prototypes
//	Functions found in Functs.cu
void arycpy(double A[],double B[],int n);
void arycpy2(double A[],double B[],int n);
void arycpy3(double A[],double B[],int n);
void set_up(double A[], double B[], double C[], int n , double h);
void get_final_q(double x[], double h, double a[], double b[], double c[], double d[], double R[], int n);
void get_final_qdot(double x[], double h, double a[], double b[], double c[], double d[], double R[], int n);
void DCAhelp(int n,bool Active_DOF[],double Coords[], double Speeds[], int dof_index[], double initZetas[],double Mass[], double Inertia[],int DOF, double Y[], int cut_off,double Body_Vectors[]);


//RK_45
//	Function used to perform Runge-Kutta 45 integration using the DCA algorithm.
//	This numerical integrator uses a fixed time step.
//		state is an array of the conditions at that timestep
//		step is the length of one timestep
//		n is the number of bodies
//		bs is a list of bodies
//		js is a list of joints
//		Y is the array where the conditions at the next timestep will be stored
//RK_45(tstep,n,Active_DOF,Coords,Speeds,dof_index,initZetas);
void RK_45(double step, int n,bool Active_DOF[],double Coords[], double Speeds[], double initZetas[],double Mass[], double Inertia[],int DOF, double Y[], double Ydot[], int dof_index[], int cut_off,double Body_Vectors[])
{
	//Comments have not yet been completed in this file because I do not know how it works
	//Variable Declarations
	double *q = new double[DOF];
	double *qdot = new double[DOF];
	double *q1 = new double[DOF];
	double *qdot1 = new double[DOF];
	double *q2 = new double[DOF];
	double *qdot2 = new double[DOF];
	double *q3 = new double[DOF];
	double *qdot3 = new double[DOF];
	double *q4 = new double[DOF];
	double *qdot4 = new double[DOF];
	double *Y1 = new double[2*n];
	double *Y2= new double[2*n];
	double *Y3=new double[2*n];
	double *Y4=new double[2*n];
	int i = 0;
	double h = step/2.0;
	
	
	arycpy(q1,Coords,DOF);
	arycpy(qdot1,Speeds,DOF);
	DCAhelp(n,Active_DOF, q1, qdot1,dof_index, initZetas,Mass,Inertia,DOF,Y1,cut_off,Body_Vectors);
	
	set_up(q,qdot1,q2,DOF,h);
	set_up(qdot,Y1,qdot2,DOF,h);
	DCAhelp(n,Active_DOF, q2, qdot2,dof_index, initZetas,Mass,Inertia,DOF,Y2,cut_off,Body_Vectors);

	set_up(q,qdot2,q3,DOF,h);
	set_up(qdot,Y2,qdot3,DOF,h);
	DCAhelp(n,Active_DOF, q3, qdot3,dof_index, initZetas,Mass,Inertia,DOF,Y3,cut_off,Body_Vectors);

	set_up(q,qdot3,q4,DOF,h*2);
	set_up(qdot,Y3,qdot4,n,h*2);
	DCAhelp(n,Active_DOF, q4, qdot4,dof_index, initZetas,Mass,Inertia,DOF,Y4,cut_off,Body_Vectors);

	get_final_q(q,h,qdot1,qdot2,qdot3,qdot4,Y, DOF);
	get_final_qdot(qdot,h,Y1,Y2,Y3,Y4,Ydot, DOF);
	
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
	delete[] Y1;
	delete[] Y2;
	delete[] Y3;
	delete[] Y4;

}

//arycpy:
//	Function used to copy array B into array A
//		A is the array where the copy will be stored
//		B is the array to be copied
//		n is the number of values in the array
void arycpy(double A[],double B[],int n)
{

	//Loop through the arrays saving every element in B into A
	for(int i =0; i<n; i++)
	{
		A[i]=B[i];
		i++;
	}
}
void set_up(double A[], double B[], double C[], int n , double h)
{
	int i = 0;	//counter

	//Loop through the elements in the arrays saving the solution in C
	while (i<n)
	{
		C[i]=A[i]+(h*B[i]);
		i++;
	}
}

//get_final_q:
//	Function used to get the final q to be used in the next timestep.
//		x is q
//		h is have of the length of one timestep
//		a is qdot1
//		b is qdot2
//		c is qdot3
//		d is qdot4
//		R is where the solution is stored
//		n is the number of bodies
void get_final_q(double x[], double h, double a[], double b[], double c[], double d[], double R[], int n)
{
	int i = 0;	//counter

	//Loop through the arrays solving for the position at the next timestep and saving it in R
	while(i<n)
	{
		R[i]= x[i] + ((h/3.0)*(a[i]+(2*b[i])+(2*c[i])+d[i]));
		i++;
	}
}

//get_final_qdot:
//	Function used to get the final qdot to be used in the next timestep.
//		x is qdot
//		h is have of the length of one timestep
//		a is qddot1
//		b is qddot2
//		c is qddot3
//		d is qddot4
//		R is where the solution is stored
//		n is the number of bodies
void get_final_qdot(double x[], double h, double a[], double b[], double c[], double d[], double R[], int n)
{
	int i = 0;	//counter

	//Loop through the arrays solving for the position at the next timestep and saving it in R
	while(i<n)
	{
		R[i]= x[i] + ((h/3.0)*(a[i]+(2*b[i])+(2*c[i])+d[i]));
		i++;
	}
}
