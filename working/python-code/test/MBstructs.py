import numpy as np
import MultiBodyFuncts as MBF

#Body is the class that defines a single body in a multibody system
class Body:
	
	#Function to create the Inertia Tensor
	def Inertia(self,I1,I2,I3):
		self.I=np.array(((I1,0,0),(0,I2,0),(0,0,I3)))
	
	#This function is not used in DCA
	def bvectors(self,qq,rr):
		self.q=qq
		self.r=rr
		self.Sq=MBF.skew(qq)
		self.Sr=MBF.skew(rr)
	
	#This function is not used in DCA
	def angvel(self,x):
		self.wn=np.dot(x,self.Pw)
		self.W0=np.transpose(MBF.skew(self.wn))
	
	#This function is not used in DCA
	def Sderiv(self):
		self.C0dot=np.dot(self.W0,self.C0)
	
	#Function to transform the Inertia Matrix with the 
	#Direction Cosine Matrices
	def Inertia2(self):
		self.I0=np.dot(np.transpose(self.C0),np.dot(self.I,self.C0))
		self.I0T=np.transpose(self.I0)
	
	#Function to create the Shifter Matrices
	def shifters(self):
		self.S01=np.zeros((6,6))
		self.S10=np.zeros((6,6))
		self.S20=np.zeros((6,6))
		self.S02=np.zeros((6,6))

		self.S02[:3,:3]=np.identity(3)
		self.S02[3:,3:]=np.identity(3)
		self.S10[:3,:3]=np.identity(3)
		self.S10[3:,3:]=np.identity(3)
		self.S20[:3,:3]=np.identity(3)
		self.S20[3:,3:]=np.identity(3)
		self.S01[:3,:3]=np.identity(3)
		self.S01[3:,3:]=np.identity(3)
				
		self.S02[:3,3:]=MBF.skew(self.r02)
		self.S10[3:,:3]=MBF.skew(self.r10)
		self.S01[:3,3:]=MBF.skew(self.r01)	
		self.S20[3:,:3]=MBF.skew(self.r20)

	#Function to Create the Mass and inverse Mass Matrices	
	def Mmatrix(self):
		self.M=np.zeros((6,6))
		self.M[:3,:3]=self.I0
		self.M[3:,3:]=np.array(((self.m,0,0),(0,self.m,0),(0,0,self.m)))		
		self.Minv=np.linalg.inv(self.M)
	
	#Function that defines position vectors between any of the
	#Three points on the body (two handles and center of mass)
	def rs(self,r,v):
		self.r10= np.dot(r*v,self.C0)
		self.r01=-1*self.r10
		self.r12= np.dot(self.l*v,self.C0)
		self.r21=-1*self.r12
		self.r20= np.dot((r)*v*-1,self.C0) 
		self.r02=-1*self.r20

	#Function that finds the initial zeta values in each timestep
	def zs(self):
		self.z11=np.dot(self.S10,np.dot(self.Minv,self.S01))
		self.z12=np.dot(self.S10,np.dot(self.Minv,self.S02))
		self.z13=np.dot(np.dot(self.S10,self.Minv),self.Fa1)
		self.z21=np.dot(self.S20,np.dot(self.Minv,self.S01))
		self.z22=np.dot(self.S20,np.dot(self.Minv,self.S02))
		self.z23=np.dot(self.S20,np.dot(self.Minv,self.Fa2))
	
	#Function to apply the state dependant forces
	#I think this is where my problem is because I am 
	#not sure if this is correct.
	def Forces(self,v):

		self.Fa1=np.zeros((6))
		self.Fa2=np.zeros((6))

		#Force from centripedal motion around the
		#body's first joint
		self.i1= np.cross(self.w,self.r10)
		self.i12=np.cross(self.w,self.r20)
		self.i2=np.cross(self.w,self.i1)
		self.i31=self.m*self.i2
		self.i32=self.m*np.cross(self.w,self.i12)
		
		
		#Total force with gravity included
		self.Fa1[3:]=9.81*self.m*v-self.i31#-self.i5
		self.Fa2[3:]=9.81*self.m*v-self.i32

#Joint is the class that defines a single joint in a multibody system		
class Joint:

	#Function to create the P and D matrices+++
	def PDmatrix(self,v,x):
		self.P=np.zeros((6,x))
		self.D=np.zeros((6,6-x))
		j=0
		k=0
		for i in range(0,6):
			if v[i]==1:
				self.P[i,j]=1
				j=j+1
			else:
				self.D[i,k] = 1
				k=k+1
			
		
