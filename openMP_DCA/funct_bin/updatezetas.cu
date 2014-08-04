#define z11 0
#define z12 6
#define z13 12
#define z21 13
#define z22 19
#define z23 25

#include <iostream>
void update(double R[6][6], double Rt[6][6], double zeta1[], double zeta2[], int n, int bodynum, int zetanum);

void MSsetup(double Minv[6][6],double S[6][6],double m, double I[3][3]);

void makeSt(double S[6][6], double r[]);

void savez3(double Minv[6][6],double S[6][6] ,double Fa[], double nZetas[], int n, int bodynum, int zetanum);

void makeR(double R[6][6], double Rt[6][6], double angle);

void printm(double A[6][6]);
void printa(double A[], int len);
void Update_Properties(double bodyZetas[],double nZetas[], int n, double state[], double m[], double l[], double II[])
{
	double R[6][6];
	double rv[3];
	double Fa[6];
	double Rt[6][6];
	double Minv[6][6];
	double S[6][6];
	double Inertia[3][3];
	double g = 9.81;	//Set the gravitational constant
	
	for(int r =0; r<6; r++)
	{
		Fa[r]=0;
	}
	
	long index1; 
	long rowlen;
	double angle;
	double w;
	angle=0;

	for(int i =0; i<n; i++)
	{

		angle=0;
		w=0;
		
		for(int c =0; c<=i; c++)
		{
			angle += state[c];
			w += state[c+n];
		}
		
		index1=i*26;
		rowlen=n*26;
		
		for(int r =0; r<3; r++)
		{
			for(int c =0; c<3; c++)
			{
				Inertia[r][c]=II[c+r*n*3+i*3];
			}
		}
		//std::cout<<(double)angle<<std::endl;
		makeR(R, Rt, angle); 
		update(R,Rt,bodyZetas,nZetas,n,i,z11);
		update(R,Rt,bodyZetas,nZetas,n,i,z12);
		update(R,Rt,bodyZetas,nZetas,n,i,z21);
		update(R,Rt,bodyZetas,nZetas,n,i,z22);
/*
		for(int r =0; r<6; r++)
		{
			for(int c =0; c<6; c++)
			{
				std::cout<<nZetas[c+r*26*n+i*26+z22]<<"  ";
			}
			std::cout<<std::endl;
		}
*/
		
		for(int r=0; r<6; r++)
	
		//Loop through the first column of S and set it to 0
		for(int r =0; r<6 ; r++)
		{
	
			nZetas[r*rowlen+index1+z13]=0;
			nZetas[r*rowlen+index1+z23]=0;
		}
		
		
		//rv is now r01
		rv[0]=-1*l[i]*sin(angle)/2;
    	rv[1]=l[i]*cos(angle)/2;
    	rv[2]=0;
				
		
		Fa[3]=-1*m[i]*w*w*rv[0];
		Fa[4]=-1*g*m[i]-m[i]*w*w*rv[1];
		
		
		MSsetup(Minv,S,m[i],Inertia);
		makeSt(S,rv);
		savez3(Minv,S,Fa,nZetas,n,i,z13);
		
		rv[0]=l[i]*sin(angle)/2;
    	rv[1]=-1*l[i]*cos(angle)/2;
    	rv[2]=0;
		Fa[3]=-1*m[i]*w*w*rv[0];
		Fa[4]=-1*g*m[i]-m[i]*w*w*rv[1];
		makeSt(S,rv);
		savez3(Minv,S,Fa,nZetas,n,i,z23);

	}
}
		










			
