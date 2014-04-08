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
		self.S0dot=np.dot(self.W0,self.S0)
	def Inertia2(self):
		self.I0=np.dot(np.transpose(self.S0),np.dot(self.I,self.S0))
		self.I0T=np.transpose(self.I0)
