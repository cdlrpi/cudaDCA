__global__ void Disassemble(float Xinv[] ,float zs[], float oldAs[] ,float newAs[], int numBlocks)
{
		const int i1 = blockIdx.x;
		const int row = threadIdx.y;
		const int col = threadIdx.x;
		//index to read Xinv[]
    	const int rXinv=col+row*numBlocks*5+i1*5;
		//indices to read zeta values
		const int r111=col+row*numBlocks*52+i1*52;
    	const int r112=r111+6;
    	const int r113=r112+6; //only col = 0 can read
    	const int r121=r113+1;
    	const int r122=r121+6;
    	const int r123=r122+6; //only col= 0 can read
    
    	const int r211=r123+1;
    	const int r212=r211+6;
    	const int r213=r212+6; //only col= 0 can read

		//index to read old acceleration and force values, only col<4 can read this
    	const int readOld=col+row*numBlocks*4+i1*4;

    	__shared__ float Vals[6][4];
    	__shared__ float temp1[6][6];
    	__shared__ float temp2[6][6];
    	__shared__ float A[6][6];

		//indices to write solution
    	const int wA11=col+row*numBlocks*8+i1*8;
    	const int wF1c1=wA11+1;
    	const int wA12=wF1c1+1;
    	const int wF1c2=wA12+1;
    	const int wA21=wF1c2+1;
    	const int wF2c1=wA21+1;
    	const int wA22=wF2c1+1;
    	const int wF2c2=wA22+1;
    	temp1[row][col]=0;
    	temp2[row][col]=0;

    	if(col<4)
    	{
    		Vals[row][col]=oldAs[readOld];
    	}
    
    	if(col==0)
    	{
    		temp1[row][col]=zs[r123];
    		temp2[row][col]=zs[r213];
   		}	
    
    	//temp1 is z23 body 1
   		//temp2 is z13 body 2
    	temp1[row][col]=temp1[row][col]-temp2[row][col];

    	//temp1 is z123-z213
    	syncthreads();
    	if(col==0)
    	{
    		temp2[row][col]=Vals[row][3];
   		}
    	//temp2 is F2c2
    	A[row][col]=zs[r212];
    	syncthreads();
    	Mat61Mult(A,temp2,temp2,row,col);
    	temp2[row][col]=-1*temp2[row][col];
    	//temp2 is -1*z212*F2c2
    	temp2[row][col]=temp2[row][col]+temp1[row][col];
    	//temp2 is -z212*F2c2+z123-z213
    	syncthreads();
    	temp1[row][col]=zs[r121];
    	//temp1 is z121
    	if(col==0)
    	{
     		A[row][col]=Vals[row][1];
    	}
    	//A is F1c1
    	syncthreads();
    	Mat61Mult(temp1,A,A,row,col);
    	//A is z121*F1c1
    	temp2[row][col]=A[row][col]+temp2[row][col];
    	//temp2 is (z121*F1c1-z212*F2c2+z123-z213)
    	if(row<5 && col<5)
    	{
    		A[row][col]=Xinv[rXinv];
    	}
    	make_W_pend(A,temp1,row,col);
    	//A is now W
    	Mat61Mult(A,temp2,A,row,col);
    	A[row][col]=-1*A[row][col];
    	//A is now -1*(W*((z121*A11)-(z212*F2c2+z123-z213)))
    	if (col==0)
    	{
    		newAs[wF1c2]=A[row][col];//F1c2
    	}
    
    	temp1[row][col]=zs[r121];
    	if(col==0)
    	{
        	temp2[row][col]=Vals[row][1];
    	}
    	syncthreads();
    	Mat61Mult(temp1,temp2,temp1,row,col);
    	temp2[row][col]=zs[r122];
    	syncthreads();
    	Mat61Mult(temp2,A,temp2,row,col);
    	temp1[row][col]=temp1[row][col]+temp2[row][col];
    	if(col==0)
    	{
        	temp2[row][col]=zs[r123];
    	}
    	syncthreads();
    	temp2[row][col]=temp1[row][col]+temp2[row][col];//temp2 should hold A2 now
    
    	if (col==0)
    	{
        	newAs[wA12]=temp2[row][col];//A12
    	}
	
    	A[row][col]=A[row][col]*-1;//Now A is the other F
    
    	if (col==0)
    	{
        	newAs[wF2c1]=A[row][col];
   		}
    
    	temp1[row][col]=zs[r211];
    	syncthreads();
    	Mat61Mult(temp1,A,A,row,col);
    	temp1[row][col]=zs[r212];
    	if(col==0)
    	{
    		temp2[row][col]=Vals[row][3];
    	}
    	syncthreads();
    	Mat61Mult(temp1,temp2,temp2,row,col);
   		if (col==0)
    	{
        	temp1[row][col]=zs[r213];
    	}
    	temp1[row][col]=A[row][col]+temp2[row][col]+temp1[row][col];//temp1 now holds A1
    	if(col==0)
    	{
        	newAs[wA21]=temp1[row][col];
        	newAs[wA11]=Vals[row][0];
        	newAs[wF1c1]=Vals[row][1];
        	newAs[wA22]=Vals[row][2];
        	newAs[wF2c2]=Vals[row][3];
    	}
}