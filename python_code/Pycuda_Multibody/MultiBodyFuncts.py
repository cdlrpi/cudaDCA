import numpy as np
#This script will hold some functions that will be used in other scripts
def directionCosine(alpha,beta,gamma):#Finds the direction cosing matrix 0S1 given by a body 1-2-3 transformation
	C=np.zeros((3,3))
	C[0,0]=np.cos(beta)*np.cos(gamma)
	C[0,1]=np.cos(alpha)*np.sin(gamma)+np.sin(alpha)*np.sin(beta)*np.cos(gamma)
	C[0,2]=np.sin(gamma)*np.sin(alpha)-np.cos(alpha)*np.sin(beta)*np.cos(gamma)
	C[1,0]=-np.cos(beta)*np.sin(gamma)
	C[1,1]= np.cos(alpha)*np.cos(gamma)-np.sin(alpha)*np.sin(beta)*np.sin(gamma)
	C[1,2]=np.sin(alpha)*np.cos(gamma)+np.cos(alpha)*np.sin(beta)*np.sin(gamma)
	C[2,0]=np.sin(beta)
	C[2,1]=-np.sin(alpha)*np.cos(beta)
	C[2,2]=np.cos(alpha)*np.cos(beta)
	return C

