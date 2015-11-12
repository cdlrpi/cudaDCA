//
//	findCutoff.cu
//
//This file contains the function that determines the most
//efficient number of assemblies to perform on the gpu

#include <cuda.h>
#include <iostream>


void cudaAssemble(double Zs[],double Xs[], int num, double nZs[], double nXs[], int odd, int newlen);
void cudaDisassemble(double OldAF[], double Zs[], double Xs[],double nZs[], double nXs[], int odd, int morelen, int lesslen, double AF[]);
void Assemble(double Zs[], double Xs[],double nZs[], double nXs[], int len, int odd, int n);
void Disassemble(double lessZs[], double lessXs[],double moreZs[], double moreXs[], double oldAs[] ,double newAs[], int num, int odd);

int findCutoff(int n, int accuracy)
{
	//Variable declarations
	int x=n;
	int numAssemblies=0;
	int count=0;
	int *numbods;
	int odd;

	//Timer creation
	float time1;
	float time2;
	time1=0;
	time2=0;
	cudaEvent_t beginEvent;
	cudaEvent_t endEvent;
	cudaEventCreate( &beginEvent );
	cudaEventCreate( &endEvent );

	//Determine the number of assemblies needed to completely assemble n bodies
	while(x!=1)
	{
		if( x%2==0)
		{
			x=x/2;
		}
		else
		{
			x++;
			x=x/2;
		}
		numAssemblies++;
	}

	//Allocate space for a matrix that holds the number of bodies at each level of assembly
	numbods=(int*)malloc(sizeof(int)*numAssemblies);
	
	//Fill numbods
	x=n;
	while(count<numAssemblies)
	{
		numbods[count]=x;
		if(x%2==0)
		{
			x=x/2;
		}
		else
		{
			x++;
			x=x/2;
		}
		count++;
	}
	count=1;
	//Begin process of finding most efficient number of assemblies
	while((count<numAssemblies) && time1>=time2)	//Compare time for gpu and cpu to complete assembly
	{
		//Create and allocate space for empty variables to mimic the dca algorithm at a level of assembly
		double* AF;
		double* Zs;
		double* Xs;
		double* nZs;
		double* nXs;
		double* AFo;
		AF=(double*)malloc(sizeof(double)*numbods[numAssemblies-1-count]*4*6);
		Zs=(double*)malloc(sizeof(double)*numbods[numAssemblies-1-count]*26*6);
		Xs=(double*)malloc(sizeof(double)*numbods[numAssemblies-1-count]*5*5);
		if(count==0)
		{
			nZs=(double*)malloc(sizeof(double)*26*6);
			nXs=(double*)malloc(sizeof(double)*25);
			AFo=(double*)malloc(sizeof(double)*6*4);
		}
		else
		{
			nZs=(double*)malloc(sizeof(double)*26*6*numbods[numAssemblies-count]);
			nXs=(double*)malloc(sizeof(double)*numbods[numAssemblies-count]*25);
			AFo=(double*)malloc(sizeof(double)*numbods[numAssemblies-count]*6*4);
		}

		//Check the parity of the number of bodies at the current level of assembly
		if(numbods[numAssemblies-1-count]%2==0)
		{
			odd=0;
		}
		else
		{
			odd=1;
		}
		
		//Check the gpu speed
		cudaEventRecord( beginEvent, 0 );	//Begin timer
		for(int i=0; i<accuracy; i++)	//Perform the operations a set number of times to vary accuracy
		{			
			//Test both assembly and disassembly at this level
			if(count==0)
			{
				cudaAssemble(Zs,Xs,numbods[numAssemblies-1-count],nZs,nXs,odd,1);
				cudaDisassemble(AFo,Zs,Xs,nZs,nXs,odd,numbods[numAssemblies-1-count],1,AF);
			}
			else
			{
				cudaAssemble(Zs,Xs,numbods[numAssemblies-1-count],nZs,nXs,odd,numbods[numAssemblies-count]);
				cudaDisassemble(AFo,Zs,Xs,nZs,nXs,odd,numbods[numAssemblies-1-count],numbods[numAssemblies-count],AF);
			}
			
		}
		//End timer
		cudaEventRecord( endEvent, 0 );
		cudaEventSynchronize( endEvent );
		cudaEventElapsedTime( &time1, beginEvent, endEvent );
		
		//Check the cpu speed
		cudaEventRecord( beginEvent,0);	//begin timing
		for(int i=0; i<accuracy; i++)	//Perfom the operations a set number of times to vary accuracy
		{
			//Test both assembly and disassembly at this level
			if(count==0)
			{
				Assemble(Zs,Xs,nZs,nXs,1,odd,numbods[numAssemblies-1-count]);
				Disassemble(nZs,nXs,Zs,Xs,AFo,AF,numbods[numAssemblies-count],odd);
			}
			else
			{
				Assemble(Zs,Xs,nZs,nXs,numbods[numAssemblies-count],odd,numbods[numAssemblies-1-count]);
				Disassemble(nZs,nXs,Zs,Xs,AFo,AF,numbods[numAssemblies-count],odd);
			}
		}
		//End timer
		cudaEventRecord( endEvent,0);
		cudaEventSynchronize(endEvent);
		cudaEventElapsedTime( &time2,beginEvent, endEvent);

		count++;
	}
	return numAssemblies-count;
}
