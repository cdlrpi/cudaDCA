import numpy as np
import MultiBodyFuncts as MBF

class Body:
	Force=np.array((0,0,0))
	def Inertia(self,I1,I2,I3):
		self.I=np.array(((I1,0,0),(0,I2,0),(0,0,I3)))

	def bvectors(self,qq,rr):
		self.q=qq
		self.r=rr
		self.Sq=MBF.skew(qq)
		self.Sr=MBF.skew(rr)

	def angvel(self,x):
		self.wn=np.dot(x,self.Pw)
		self.W0=np.transpose(MBF.skew(self.wn))

	def Sderiv(self):
		self.C0dot=np.dot(self.W0,self.C0)
	def Inertia2(self):
		self.I0=np.dot(np.transpose(self.C0),np.dot(self.I,self.C0))
	
		self.I0T=np.transpose(self.I0)
	def shifters(self):
		self.S12=np.zeros((6,6))
		self.S10=np.zeros((6,6))
		self.S20=np.zeros((6,6))
		self.S21=np.zeros((6,6))

		self.S12[:3,:3]=np.identity(3)
		self.S12[3:,3:]=np.identity(3)
		self.S10[:3,:3]=np.identity(3)
		self.S10[3:,3:]=np.identity(3)
		self.S20[:3,:3]=np.identity(3)
		self.S20[3:,3:]=np.identity(3)
		self.S21[:3,:3]=np.identity(3)
		self.S21[3:,3:]=np.identity(3)
		self.S12[:3,3:]=MBF.skew(self.r12)
		self.S10[:3,3:]=MBF.skew(self.r10)
		self.S21[:3,3:]=MBF.skew(self.r21)	
		self.S20[:3,3:]=MBF.skew(self.r20)

	def Mmatrix(self):
		self.M=np.zeros((6,6))
		self.M[0:3,0:3]=self.I0
		self.M[3:,:3]=np.identity(3)
		self.M[:3,3:]=np.identity(3) #this is how i've been forcing it to be invertible
		self.M[3:,3:]=np.array(((self.m,0,0),(0,self.m,0),(0,0,self.m)))
	
		self.Minv=np.zeros((6,6))
		self.Minv=np.linalg.inv(self.M)
	def rs(self,r,v):
		self.r10= np.dot(r*v,self.C0)
		self.r12= np.dot(self.l*v,self.C0)
		self.r21=-1*self.r12
		self.r20= np.dot((self.l-r)*v*-1,self.C0) 
	def zs(self):
		self.z11=self.Minv
		self.z12=np.dot(self.Minv,self.S12)
		self.z13=np.dot(self.Minv,np.dot(self.S10,self.Fa))
		self.z21=np.dot(self.Minv,self.S21)
		self.z22=self.Minv
		self.z23=np.dot(self.Minv,np.dot(self.S20,self.Fa))
	
	def gravity(self,v):
		self.Fa=np.zeros((6))
		self.Fa[3:]=9.81*self.m*v



class Joint:
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
			
		
