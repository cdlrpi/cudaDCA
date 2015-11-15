//	cudaDCA.cu
//
//This file contains the recursive DCA function, and the function that is used to invoke DCA and
//interperate the results.

//Included Files

#include <iostream>

//Function Prototypes
//	Functions found in this file
void RecDCA(double Zs[], int n, int i, double AF[], int cut_off,double Xs[],float times[]);

//	Functions found in Init_setup.cu
void CudaInitialize(double m[], double l[], double I[], double x[], int n, double Zs[]);

//	Functions found in Assemble_setup.cu
void cudaAssemble(double Zs[],double Xs[], int num, double nZs[], double nXs[], int odd, int newlen,float times[]);

//	Functions found in Disassemble_setup.cu
void cudaDisassemble(double OldAF[], double Zs[], double Xs[],double nZs[], double nXs[], int odd, int morelen, int lesslen, double AF[],float times[]);

//	Functions found in Assemble.cu
void Assemble(double Zs[], double Xs[],double nZs[], double nXs[], int len, int odd, int n);

//	Functions found in Disassemble.cu
void Disassemble(double lessZs[], double lessXs[],double moreZs[], double moreXs[], double oldAs[] ,double newAs[], int num, int odd);

//	Functions found in SolveBCs.cu
void solve_BCs(double Zs[], double Xs[], double AF[]);

void printa(double A[], int n);
void printm(double A[6][6]);

void Update_Properties(double bodyZetas[],double nZetas[], int n, double state[], double m[], double l[], double II[]);
//DCAhelp:
//	Function that prepares the list of bodies for DCA and finds the final state vector
//		state is the state of the system at that timestep
//		bs is a list of bodies used for initialization
//		js is a list of joints
//		n is the number of bodies
//		Y is the array where the final velocities and accelerations are stored
void DCAhelp(double state[], double m[], double l[], double I[],int n, double Y[],int cut_off, double bZs[], float times[], int reps)
{
	//Create the list that will hold all acceleration and force values for all bodies
	double *AF = (double*) malloc(sizeof(double)*n*4*6);
	//double A[6][n*2];	//Create the matrix where only the accelerations will be stored
	double *Zs=(double*)malloc(sizeof(double)*n*6*26);
	double *Xs=(double*)malloc(sizeof(double)*n*5*5);
	

	Update_Properties(bZs,Zs,n,state,m,l,I);

	//CudaInitialize(m,l,I, state, n, Zs);	//Initialize the bodies, finding all zeta values
	
	//Pass the list of bodies to DCA and return the accelerations 
	//and forces of both handles of every body in the list

	RecDCA(Zs, n, 0, AF, cut_off,Xs,times);

	Y[n]=AF[8*n];	//For a pendulum, the fist acceleration value is in A[2][0]

	for(int i = n+1, j=2; i<n*2; i++, j+=2)	//Loop through the acceleration matrix
	{
		Y[i]=AF[2*4*n+2*j]-AF[2*4*n+2*(j-1)]; //Find and save all generalized accelerations
	}

	for(int i = 0; i<n; i++)	//Loop through the state vector
	{
		Y[i]=state[i+n];	//Save the velocities
	}

	//Free memory
	free(AF);
	free(Zs);
	free(Xs);		
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
void RecDCA(double Zs[], int n, int i, double AF[], int cut_off,double Xs[], float times[])
{
	if (n==1)	//If there is only 1 body
	{	
		
		solve_BCs(Zs,Xs, AF);	//Solve the boundary conditions and find the acceleratins and forces
	
		
	}
	else	//If there is more than 1 body
	{	
		int newlen;	//New number of bodies after assembly
		int odd = 0;	//Flag to keep track of the parity of the length of the list of
		
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
		double *nZs=(double*)malloc(sizeof(double)*newlen*26*6);
		double *nXs=(double*) malloc(sizeof(double)*(newlen)*5*5);
		if(n>cut_off)
		{
			//std::cout<<'g';
			//std::cout<<n-cut_off;
			//std::cout<<"YYY"<<cut_off<<"YYY";
			//std::cout<<"III"<<n<<"III";
			cudaAssemble(Zs,Xs, n, nZs, nXs , odd, newlen, times);	//Assemble the bodies, storing them in newbds		
			
		}
		else
		{
		//	std::cout<<"cpu"<<std::endl;
		
 	
			Assemble(Zs,Xs,nZs,nXs, newlen,odd, n);	//Assemble the bodies, storing them in newbds
		}
		//Create a list of accelerations and forces of the new bodies.
		double *AFo=(double*)malloc(sizeof(double)*6*newlen*4);

		//Call the DCA function again to return the accelerations and forces of the new bodies
		RecDCA(nZs,newlen,i+1 ,AFo,cut_off,nXs,times);

		//Knowing the accelerations and forces of the new bodies, the new bodies can be disassembled
		//again, finding the accelerations and forces of the old bodies.
		if(n>cut_off)
		{
			
			cudaDisassemble(AFo, Zs,Xs , nZs,nXs, odd, n, newlen, AF, times);
			
		}
		else
		{
			Disassemble(nZs,nXs,Zs,Xs,AFo, AF, newlen,odd);
		}
		//Free memory
		free(nZs);
		
		free(nXs);
		
		free(AFo);
		 
	
	}
}
		
			


