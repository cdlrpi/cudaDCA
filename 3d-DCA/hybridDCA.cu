//	cudaDCA.cu
//
//This file contains the recursive DCA function, and the function that is used to invoke DCA and
//interperate the results.

//Included Files

#include <iostream>

//Function Prototypes
//	Functions found in this file
void RecDCA(double Zs[], int n, int i, double AF[], int cut_off,double Xs[]);

//	Functions found in Init_setup.cu
void CudaInitialize(double m[], double l[], double I[], double x[], int n, double Zs[]);

//	Functions found in Assemble_setup.cu
void cudaAssemble(double Zs[],double Xs[], int num, double nZs[], double nXs[], int odd, int newlen);

//	Functions found in Disassemble_setup.cu
void cudaDisassemble(double OldAF[], double Zs[], double Xs[],double nZs[], double nXs[], int odd, int morelen, int lesslen, double AF[]);

//	Functions found in Assemble.cu

//	Functions found in SolveBCs.cu
void solve_BCs(double Zs[], double Xs[], double AF[]);

void printa(double A[], int n);
void printm(double A[6][6]);
void update(double Mass[], double Inertia[], double Init_Zetas[], double Zetas[], int n, double Body_Vectors[], double ang_vels[], double DCMs[],double Ps[], double PdotUs[],double Ds[],bool active_DOF[], double speeds[],int dof_index[]);
void kinematics(bool Active_DOF[], double Coords[], double DCMs[], double speeds[], double omegas[], int dof_index[],int n);
void getGenAccel(double Y[], double AF[],double omegas[], bool Active_DOF[], int n,double PdotUs[], double Ps[]);
void multPinv(double Ps[], double tempy[], int b, int n);
void RecDCA(double Zs[], int n, int i, double AF[], int cut_off,double Xs[], double DCMs[], double omegas[], int dof_index[], bool Active_DOF[], double Speeds[],double PdotUs[], double Ds[], int Pindex[],int numbods);
void Assemble(double Zs[], double Xs[],double nZs[], double nXs[], int len, int odd, int n, double PdotUs[], double Ds[], int Pindex[],int numbods, bool active_DOF[]);
void Disassemble(double lessZs[], double lessXs[],double moreZs[], double moreXs[], double oldAs[] ,double newAs[], int num, int odd,double Ds[], int Pindex[],int numbods);
void multPinv(double Ps[], double tempy[], int b, int n);
void Mat61Mult2(double A[6][6], double B[6], double C[6]);
//DCAhelp:
//	Function that prepares the list of bodies for DCA and finds the final state vector
//		state is the state of the system at that timestep
//		bs is a list of bodies used for initialization
//		js is a list of joints
//		n is the number of bodies
//		Y is the array where the final velocities and accelerations are stored
void DCAhelp(int n,bool Active_DOF[],double Coords[], double Speeds[], int dof_index[], double initZetas[],double Mass[], double Inertia[],int DOF, double Y[], int cut_off, double Body_Vectors[])
{
	double *Zetas = new double[n*26*6];
	double *Xs = new double[n*5*5];
	double *DCMs = new double[n*3*3];
	double *omegas= new double[n*3];
	double *AF = new double[n*4*6];
	double *Ps = new double[n*6*6];
	double *PdotUs = new double[n*6];
	double *Ds = new double[n*6*6];
	int *Pindex = new int[n];
	for(int i =0; i<n; i++)
	{
		Pindex[i]=i;
	}

	kinematics(Active_DOF, Coords, DCMs, Speeds, omegas, dof_index, n);
	update(Mass, Inertia, initZetas, Zetas, n, Body_Vectors, omegas, DCMs,Ps,PdotUs,Ds,Active_DOF,Speeds,dof_index);
	RecDCA(Zetas, n, 0, AF, cut_off,Xs,DCMs,omegas,dof_index,Active_DOF,Speeds,PdotUs,Ds,Pindex,n);
	getGenAccel(Y,AF,omegas,Active_DOF,n,PdotUs,Ps);

	//Free memory
	delete[] Zetas;
	delete[] Xs;
	delete[] DCMs;
	delete[] omegas;
	delete[] AF;
	delete[] Ps;
	delete[] PdotUs;
	delete[] Ds;
	delete[] Pindex;
}


void getGenAccel(double Y[], double AF[],double omegas[], bool Active_DOF[], int n, double PdotUs[], double Ps[])
{
	int c =0;
	//double A1[6];
	//double A2[6];
	double tempy[6];

	for(int i =0; i<6; i++)
	{
		tempy[i]=AF[i*4]-PdotUs[i];
	}
	multPinv(Ps,tempy,0,n);
	for(int i =0; i<6; i++)
	{
		if(Active_DOF[i])
		{
			Y[c]=tempy[i];
			c++;
		}
	}	
	for(int i =1; i<n; i++)
	{
		for(int k =0; k<6; k++)
		{
			tempy[k]=AF[i*4+k*4*n+1]-AF[i*4+k*4*n-1]-PdotUs[i*6+k];
		}
		multPinv(Ps,tempy,i,n);
		for(int j=0; j<6; j++)
		{
			if(Active_DOF[6*i+j])
			{
				Y[c]=tempy[j];
				c++;
			}
		}
	}
}

void multPinv(double Ps[], double tempy[], int b, int n)
{
	double Pinv[6][6];
	
	for( int r =0; r<6; r++)
	{
		for(int c =0; c<6; c++)
		{
			Pinv[c][r]=Ps[c+b*6+r*6*n];
		}
	}
	Mat61Mult2(Pinv,tempy,tempy);
}
	
			
//RecDCA:
//	Function used to solve for the velocty and acceleration of the list of bodies at
//	the current timestep.  This is a recursive function that continues to call itself
//	until there is a single body left.  Once at this point the accelerations and forces
//	are found using the boundary conditions of a pendulum.  These values are then returned
//	to the previous level of recursion which then finds the new accelerations and forces
//	for the disassembled bodies.  This continues until all bodies are disassembled, ultimately
//	returning the forces and accelerations at both handles of every body in the system.  These
//	results are intererated by DCAhelp (above) to obtain the actual generalized accelerations.  
//		bodies is the list of bodies
//		n is the number of bodies
//		i is the level of recursion
//		AF is the array in which the accelerations and forces at the handles of the bodies 
//			will be stored.
void RecDCA(double Zs[], int n, int i, double AF[], int cut_off,double Xs[], double DCMs[], double omegas[], int dof_index[], bool Active_DOF[], double Speeds[],double PdotUs[], double Ds[], int Pindex[],int numbods)
{
	if (n==1)	//If there is only 1 body
	{	
		
		solve_BCs(Zs,Xs, AF);	//Solve the boundary conditions and find the acceleratins and forces
	
		
	}
	else	//If there is more than 1 body
	{	
		int newlen;	//New number of bodies after assembly
		bool odd = 0;	//Flag to keep track of the parity of the length of the list of
		
		if(n % 2 == 0)	//If there is an even number of bodies
		{
			newlen = (int) (n/2);	//The new number of bodies will be half the original number
		}
		else	//If there is an odd number of bodies
		{
			newlen = (int)((n+1)/2);	//The new number of bodies will be half the original number
										//rounded down, plus 1
			odd = 1;	//odd is set to 1 because there are an odd number of bodies
		}
		double *nZs=new double[newlen*26*6];
		double *nXs=new double[newlen*5*5];
		int *nPindex = new int[newlen];

		for(int k =0; k<newlen-odd; k++)
		{
			nPindex[i]=Pindex[i*2];
		}
		if(odd)
		{
			nPindex[newlen-1]=Pindex[n-1];
		}

/*
		if(i<cut_off)
		{
			cudaAssemble(Zs,Xs, n, nZs, nXs , odd, newlen);	//Assemble the bodies, storing them in newbds			
		}
		else
		{*/
 	
			Assemble(Zs,Xs,nZs,nXs, newlen,odd, n,PdotUs,Ds,Pindex,numbods,Active_DOF);	//Assemble the bodies, storing them in newbds
		//}

		//Create a list of accelerations and forces of the new bodies.
		double *AFo= new double[6*newlen*4];

		//Call the DCA function again to return the accelerations and forces of the new bodies
		RecDCA(nZs,newlen,i+1 ,AFo,cut_off,nXs,DCMs,omegas,dof_index,Active_DOF,Speeds,PdotUs,Ds,nPindex,numbods);

		//Knowing the accelerations and forces of the new bodies, the new bodies can be disassembled
		//again, finding the accelerations and forces of the old bodies.
		/*if(i<cut_off)
		{
			
			cudaDisassemble(AFo, Zs,Xs , nZs,nXs, odd, n, newlen, AF);
			
		}
		else
		{*/
		Disassemble(nZs,nXs,Zs,Xs,AFo, AF, newlen,odd,Ds,Pindex,numbods);
		//}
		//Free memory
		delete[] nZs;
		
		delete[] nXs;
		
		delete[] AFo;

		delete[] nPindex;
		 
	
	}
}
		
			

