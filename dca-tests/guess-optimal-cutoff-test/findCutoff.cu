#include "deviceAssemble.h"
#include "deviceDisassemble.h"
#include <cuda.h>
#include <iostream>
#include "deviceFuncts.h"

void cudaAssemble(double Zs[],double Xs[], int num, double nZs[], double nXs[], int odd, int newlen);
void cudaDisassemble(double OldAF[], double Zs[], double Xs[],double nZs[], double nXs[], int odd, int morelen, int lesslen, double AF[]);
void Assemble(double Zs[], double Xs[],double nZs[], double nXs[], int len, int odd, int n);
void Disassemble(double lessZs[], double lessXs[],double moreZs[], double moreXs[], double oldAs[] ,double newAs[], int num, int odd);

int main(void)
{

	int n=3200;
	int accuracy = 20;
	int x=n;
	int numAssemblies=0;
	int count=0;
	int *numbods;
	int odd;
	float time1;
	float time2;
	time1=0;
	time2=0;
	cudaEvent_t beginEvent;
	cudaEvent_t endEvent;
	cudaEventCreate( &beginEvent );
	cudaEventCreate( &endEvent );

	
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
	numbods=(int*)malloc(sizeof(int)*numAssemblies);
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
	for(int i=0; i<numAssemblies; i++)
	{
		std::cout<<numbods[i]<<"  ";
	}
	std::cout<<std::endl;
	//std::cout<<"here"<<std::endl;
	while((count<numAssemblies) && time1>=time2)
	{
		//std::cout<<"here"<<std::endl;
		//std::cout<<numbods[numAssemblies-1-count];
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
		if(numbods[numAssemblies-1-count]%2==0)
		{
			odd=0;
		}
		else
		{
			odd=1;
		}
		cudaEventRecord( beginEvent, 0 );

		for(int i=0; i<accuracy; i++)
		{			
			//Test Assembly
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
		cudaEventRecord( endEvent, 0 );
		cudaEventSynchronize( endEvent );
		cudaEventElapsedTime( &time1, beginEvent, endEvent );
		
		cudaEventRecord( beginEvent,0);
	
		for(int i=0; i<accuracy; i++)
		{
			
			if(count==0)
			{
				//std::cout<<"k"<<std::endl;
				Assemble(Zs,Xs,nZs,nXs,1,odd,numbods[numAssemblies-1-count]);
				//std::cout<<"A"<<std::endl;
				Disassemble(nZs,nXs,Zs,Xs,AFo,AF,numbods[numAssemblies-count],odd);
			}
			else
			{
				Assemble(Zs,Xs,nZs,nXs,numbods[numAssemblies-count],odd,numbods[numAssemblies-1-count]);
				
				Disassemble(nZs,nXs,Zs,Xs,AFo,AF,numbods[numAssemblies-count],odd);
			}
		}
		
		cudaEventRecord( endEvent,0);
		cudaEventSynchronize(endEvent);
		cudaEventElapsedTime( &time2,beginEvent, endEvent);
		//std::cout<<time1<<std::endl;
		count++;
	}
	//std::cout<<"okay"<<std::endl;
	std::cout<<"The best number of gpu assemblies is:\t" << numAssemblies-count<<std::endl;
	return 0;
}
