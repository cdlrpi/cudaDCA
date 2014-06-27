//	deviceAssemble.h
//
//This file contains the function that assembles bodies

//Function Prototypes
//	Functions found in deviceFuncts.h
__device__ void make_W(float Xinv[6][6], float D[6][6],int row , int col);                    
__device__ void get_X(float z1[6][6], float z2 [6][6], float D[6][6],int row , int col);
__device__ void invert_X(float X[6][6], float A[6][6],int row, int col);
__device__ void Mat61Mult(float A[6][6], float B[6][6], float C[6][6],int row, int col);
__device__ void Mat66Mult(float A[6][6], float B[6][6], float C[6][6],int row, int col);

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
__global__ void Assemble(float oldZs[],float newZs[],float Xinvs[], int nn)
{
	//Three shared memory matrices used for matrix operations	
	__shared__ float z1[6][6];
	__shared__ float z2[6][6];
	__shared__ float A[6][6];

	//Variables to distinguish between threads
	const int i1 = blockIdx.x;
	const int row = threadIdx.y;
	const int col = threadIdx.x;

	//Indices used to read the zetas of the bodies being assembled
	const int r111=col+row*nn*52+i1*52;
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
	const int i11 = col+row*rowlen+i1*26;
	const int i12=i11+6;
	const int i13=i12+6; //only col=0 can write
	const int i21=i13+1;
	const int i22=i21+6;
	const int i23=i22+6; //only col = 0 can write

	//index to save X inverse matrix into Xinv
	const int wXinv = col+row*(nn)*5+i1*5;//only write if col !=5 and row !=5
	__syncthreads();

   	//z11 of body 2 and z22 of body 1 are obtained to construct X
	z1[row][col]=oldZs[r211];	//z11 of body 2
	z2[row][col]=oldZs[r122];	//z22 of body 1
	__syncthreads();//To ensure matrices are done loading

	get_X(z1,z2,A,row,col);//Get the X matrix and save it in A
	invert_X(A,z1,row,col);//Invert X and save it in A

	if (col!=5 && row !=5)//Xinv is only a 5x5
	{
    	Xinvs[wXinv]=A[row][col];//Save Xinv
	}
	make_W(A,z1,row,col);//Make W out of X and save it in A
	//A now holds W

	z1[row][col]=oldZs[r112];//Load z12 for body 1 into z1
	z2[row][col]=oldZs[r212];//Load z12 for body 2 into z2
	__syncthreads();//To ensure matrices are done loading

	Mat66Mult(A,z2,z2,row,col);
	Mat66Mult(z1,z2,z2,row,col);//z2 now holds the new z12
	newZs[i12]=z2[row][col];//Save the new z12

	z2[row][col]=oldZs[r121];//Load z21 for body 1 into z2
	__syncthreads();//To ensure the matrix is done loading

	Mat66Mult(A,z2,z2,row,col);
	Mat66Mult(z1,z2,z2,row,col);
	z1[row][col]=oldZs[r111];//z1 now holds z11 for body 1
	__syncthreads();//To ensure the matrix is done loading

	z2[row][col]=z1[row][col]-z2[row][col];//z2 now holds the new z11
	newZs[i11]=z2[row][col];//Save the new z11

	z1[row][col]=oldZs[r221];//z21 for body 2
	z2[row][col]=oldZs[r121];//z21 for body 1
		__syncthreads();//To ensure matrices are done loading

	Mat66Mult(A,z2,z2,row,col);
	Mat66Mult(z1,z2,z2,row,col);//z2 now holds the new z21
	newZs[i21]=z2[row][col];//Save the new z21

	z2[row][col]=oldZs[r212];//z12 for body 2
	__syncthreads();//To ensure the matrix is done loading

	Mat66Mult(A,z2,z2,row,col);
	Mat66Mult(z1,z2,z2,row,col);
	z1[row][col]=oldZs[r222];//z22 for body 2
	__syncthreads();//To ensure the matrix is done loading

	z2[row][col]=z1[row][col]-z2[row][col];//z2 now holds the new z22
	newZs[i22]=z2[row][col];//Save the new z22

	//now find Y
	if(col==0)
	{
		z1[row][col]=oldZs[r123];//z23 for body 1
		z2[row][col]=oldZs[r213];//z13 for body 2
	}
	__syncthreads();//To ensure the matrices are done loading

	z1[row][col]=z1[row][col]-z2[row][col];
	Mat61Mult(A,z1,A,row,col);//col = 0 of A now holds Y

	z2[row][col]=oldZs[r112];//z12 for body 1
	if (col==0)
	{
		z1[row][col]=oldZs[r113];//z13 for body 1
	}
	__syncthreads();//To ensure the matrix is done loading

	Mat61Mult(z2,A,z2,row,col);
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
	__syncthreads();//To ensure the matrices are done loading

	Mat61Mult(z2,A,z2,row,col);
	z1[row][col]=z1[row][col]+z2[row][col];//z1 now holds the new z23
	if(col==0)
	{
	    newZs[i23]=z1[row][col];//Save the new z23
	}
}
