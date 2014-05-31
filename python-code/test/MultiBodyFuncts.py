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

def simdcm(angle,vector):
	C=np.zeros((3,3))
	
	cos=np.cos(angle)
	sin=np.sin(angle)
	C[0,0]=C[1,1]=C[2,2]=cos
	C[0,1]=C[0,2]=C[1,2]=C[2,1]=C[2,0]=C[1,0]=sin
	if vector[0] == 1:
		C[0,0]=1
		C[0,2]=C[0,1]=C[1,0]=C[2,0]=0
		C[2,1]*=-1
		return C
	if vector[1]==1:
		C[1,1]=1
		C[0,1]=C[1,0]=C[1,2]=C[2,1]=0
		C[0,2]*=-1
		return C
	else:
		C[2,2]=1
		C[0,2]=C[1,2]=C[2,0]=C[2,1]=0
		C[1,0]*=-1
		return C


def skew(V):
	one=V[0]
	two=V[1]
	three=V[2]
	return np.array([[0,-1*three,two],[three,0,-1*one],[-1*two,one,0]])

def vb2pts(va,awb,rab):
	return (va+np.cross(awb,rab))

def ab2pts(aa,aab,rab,awb):
	return aa+np.cross(aab,rab)+np.cross(awb,np.cross(awb,rab))
def v1pt(vbar,v):
	return vbar+v
def a1pt(abar,a,awb,v):
	return abar+a+2*np.cross(awb,v)
def genvel(v1,v2,awb,rab):
	return v1+v2+np.cross(awb,rab)
def genacc(a1,a2,aab,rab,awb, v2):
	return a1+a2+np.cross(aab,rab)+np.cross(awb,np.cross(awb,rab))+2*np.cross(awb,v2)

def PendEnergy(Y,bodies):
	val1=0
	val2=0
	height1=0
	temp_V=0
	height2=0
	n=int(len(Y[0,:])/2)
	theta=0
	w=0
	r=np.zeros((2))
	V_cm=np.zeros((2))
	h=0
	PE=np.zeros((len(Y[:,0])))
	KEr=np.zeros_like(PE)
	KEt=np.zeros_like(PE)
	energy=np.zeros((len(Y[:,0]),3))
	for i in range (0,len(Y[:,0])):
		val1=0
		val2=0
		height1=0
		temp_V=0
		height2=0
		theta=0
		w=0
		h=0
		for j in range (0,n):
			val1+=Y[i,j]
			val2+=Y[i,j+n]
			theta=val1
			w=val2
			I=bodies[j].m*(bodies[j].l**2)/12
			if j == 0:
				r[:]=np.array((bodies[j].l*np.cos(theta)/2,bodies[j].l*np.sin(theta)/2))
				V_cm[:]=r[:]*w
				temp_V+=w*r*2
				h=(bodies[0].l/2)-r[0]
				height1+=bodies[0].l
				height2+=r[0]*2
			else:

				r[:]=np.array((bodies[j].l*np.cos(theta)/2,bodies[j].l*np.sin(theta)/2))
				V_cm[:]=temp_V+w*r[:]
				temp_V+=w*r*2
				h=height1+(bodies[j].l/2)-r[0]-height2
				height1+=bodies[j].l
				height2+=r[0]*2
			PE[i]+=bodies[j].m*9.81*h
			KEr[i]+=.5*I*w*w
			KEt[i]+=.5*bodies[j].m*np.dot(V_cm,V_cm)
			
	KE=KEr+KEt
	TE=PE+KE
	energy[:,0]=KE
	energy[:,1]=PE
	energy[:,2]=TE
	return energy
				
