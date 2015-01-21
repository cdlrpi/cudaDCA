//	deviceAssemble.h
//
//This file contains the function that assembles bodies

//Function Prototypes
//	Functions found in deviceFuncts.h
__device__ void make_W_gpu(double Xinv[6][6], double D[6][6],int row);                    
__device__ void get_X_gpu(double z1[6][6], double z2 [6][6], double D[6][6],int row);
__device__ void invert_X_gpu(double X[6][6], double A[6][6],int row);
__device__ void Mat61Mult_gpu(double A[6][6], double B[6][6], double C[6][6],int row);
__device__ void Mat66Mult_gpu(double A[6][6], double B[6][6], double C[6][6],int row);
//__device__ void printit_gpu(double A[6][6], int n);
//Assemble:
//	Function used to assemble a list of bodies into a list of bodies that is 
//	half the size of the original list. This function assumes that the number of 
//	bodies it must assemble is an even number.  The list of old zeta values 
//	is cycled through looking at two bodies at a time.  These two bodies are assembled
//	into one new body, saving the zeta values for the new body in newZs.
//		oldZs is the list of zeta values for the bodies to be assembled
//		newZs is the list of zeta values for the newly made bodies
//		Xinvs is the list of X inverse matrices for teh newly made bodies
//		nn is the number of new bodies
__global__ void Assemble_gpu(double oldZs[],double newZs[],double Xinvs[], int nn, int olen)
{

	//Three shared memory matrices used for matrix operations	
	__shared__ double z1[6][6];
	__shared__ double z2[6][6];
	__shared__ double A[6][6];

	//Variables to distinguish between threads
	const int i1 = blockIdx.x;
	const int row = threadIdx.y;
	//const int col = threadIdx.x;

	//Indices used to read the zetas of the bodies being assembled
	const int r111=row*olen*26+i1*52;
	const int r112=r111+6;
	const int r113=r112+6; //only col = 0 can read
	const int r121=r113+1;
	const int r122=r121+6;
	const int r123=r122+6; //only col= 0 can read

	const int r211=r123+1;
	const int r212=r211+6;
	const int r213=r212+6; //only col= 0 can read
	const int r221=r213+1; 
	const int r222=r221+6;
	const int r223=r222+6; //only col= 0 can read

	
	int rowlen;	//Length of one row of the the new zeta matrix
	rowlen = (nn)*26;

	//indices to write the zeta values of the newly assembled body
	const int i11 = row*rowlen+i1*26;
	const int i12=i11+6;
	const int i13=i12+6; //only col=0 can write
	const int i21=i13+1;
	const int i22=i21+6;
	const int i23=i22+6; //only col = 0 can write

	//index to save X inverse matrix into Xinv
	const int wXinv = row*(nn)*5+i1*5;//only write if col !=5 and row !=5
	
	__syncthreads();

   	//z11 of body 2 and z22 of body 1 are obtained to construct X
	for(int col = 0; col<6; col++)
	{
		z1[row][col]=oldZs[r211+col];	//z11 of body 2
		z2[row][col]=oldZs[r122+col];	//z22 of body 1
	}
	__syncthreads();//To ensure matrices are done loading
	
	get_X_gpu(z1,z2,A,row);//Get the X matrix and save it in A
	invert_X_gpu(A,z1,row);//Invert X and save it in A
	//printit_gpu(A, 13);
	//printit_gpu(z2, 13);
	if (row !=5)//Xinv is only a 5x5
	{
		for(int col=0; col<5;col++)
		{
    		Xinvs[wXinv+col]=A[row][col];//Save Xinv
		}
	}
	make_W_gpu(A,z1,row);//Make W out of X and save it in A
	//A now holds W
	for(int col=0; col<6; col++)
	{
		z1[row][col]=oldZs[r112+col];//Load z12 for body 1 into z1
		z2[row][col]=oldZs[r212+col];//Load z12 for body 2 into z2
	}
	__syncthreads();//To ensure matrices are done loading

	Mat66Mult_gpu(A,z2,z2,row);
	Mat66Mult_gpu(z1,z2,z2,row);//z2 now holds the new z12
	for(int col=0; col<6; col++)
	{
		newZs[i12+col]=z2[row][col];//Save the new z12
		z2[row][col]=oldZs[r121+col];//Load z21 for body 1 into z2
	}
	__syncthreads();//To ensure the matrix is done loading

	Mat66Mult_gpu(A,z2,z2,row);
	Mat66Mult_gpu(z1,z2,z2,row);
	for(int col=0; col<6; col++)
	{
		z1[row][col]=oldZs[r111+col];//z1 now holds z11 for body 1
	}
	__syncthreads();//To ensure the matrix is done loading

	for(int col=0; col<6; col++)
	{
		z2[row][col]=z1[row][col]-z2[row][col];//z2 now holds the new z11
		newZs[i11+col]=z2[row][col];//Save the new z11

		z1[row][col]=oldZs[r221+col];//z21 for body 2
		z2[row][col]=oldZs[r121+col];//z21 for body 1
	}
		__syncthreads();//To ensure matrices are done loading

	Mat66Mult_gpu(A,z2,z2,row);
	Mat66Mult_gpu(z1,z2,z2,row);//z2 now holds the new z21
	for(int col=0; col<6; col++)
	{
		newZs[i21+col]=z2[row][col];//Save the new z21

		z2[row][col]=oldZs[r212+col];//z12 for body 2
	}
	__syncthreads();//To ensure the matrix is done loading

	Mat66Mult_gpu(A,z2,z2,row);
	Mat66Mult_gpu(z1,z2,z2,row);
	for(int col=0; col<6; col++)
	{
		z1[row][col]=oldZs[r222];//z22 for body 2
	}
	__syncthreads();//To ensure the matrix is done loading
	for(int col=0; col<6; col++)
	{
		z2[row][col]=z1[row][col]-z2[row][col];//z2 now holds the new z22
		newZs[i22+col]=z2[row][col];//Save the new z22
	}
	//now find Y
	
		z1[row][0]=oldZs[r123];//z23 for body 1
		z2[row][0]=oldZs[r213];//z13 for body 2
	
	__syncthreads();//To ensure the matrices are done loading

	for(int col=0; col<6; col++)
	{
		z1[row][col]=z1[row][col]-z2[row][col];
	}
	Mat61Mult_gpu(A,z1,A,row);//col = 0 of A now holds Y

	for(int col=0; col<6; col++)
	{
		z2[row][col]=oldZs[r112];//z12 for body 1
		if (col==0)
		{
			z1[row][col]=oldZs[r113];//z13 for body 1
		}
	}
	__syncthreads();//To ensure the matrix is done loading

	Mat61Mult_gpu(z2,A,z2,row);
	for(int col=0; col<6; col++)
	{
		z1[row][col]=z1[row][col]-z2[row][col];//z1 now holds the new z13
		if(col==0)
		{
			newZs[i13]=z1[row][col];//Save the new z13
		}
		z2[row][col]=oldZs[r221];//z21 for body 2
		if(col==0)
		{
			z1[row][col]=oldZs[r223];//z23 for body 2
		}
	}
	__syncthreads();//To ensure the matrices are done loading

	Mat61Mult_gpu(z2,A,z2,row);
	for(int col=0; col<6; col++)
	{
	z1[row][col]=z1[row][col]+z2[row][col];//z1 now holds the new z23
		if(col==0)
		{
			newZs[i23]=z1[row][col];//Save the new z23
		}
	}


}
