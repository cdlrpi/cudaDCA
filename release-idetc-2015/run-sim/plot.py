

import numpy as np
import math
import matplotlib.pyplot as plt
class Body:
	pass



#Y = np.load('mat.npy')
with open('Solution.mtx') as file:
    YY = [[float(digit) for digit in line.split()] for line in file]

with open('Body_Values.mtx') as file:
    vVals = [[float(digit) for digit in line.split()] for line in file]

Vals=np.array(vVals,dtype=np.float64)
Y=np.array(YY,dtype=np.float64)
val1=0
val2=0
height1=0
temp_V=0
height2=0
n=int(len(Y[0,:])/2)
bodies=[]
for i in range(0,n):
	bodies.append(Body())
	bodies[i].m=Vals[0,i+1];
	bodies[i].l=Vals[1,i+1];
Time = np.linspace(0,Vals[1,0],len(Y[:,0]))
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
print(len(Time)-len(TE))
plt.plot(Time,TE-TE[0])
plt.xlabel("Time [s]")
plt.ylabel("energy")
plt.title("System Energy")
plt.show()

plt.plot(Time,PE,Time,KE)
plt.xlabel("Time[s]")
plt.ylabel("energy")
plt.title("Kinetic and Potential Energy")
plt.show()

# Solution Plot
#--------------------
plt.plot(Time,Y[:,:n])
plt.xlabel("Time [s]")
plt.ylabel("Generalized Coordinates [Rad]")
plt.title("System Response")

plt.show()

plt.plot(Time,Y[:,n:])
plt.xlabel(("Time[s]"))
plt.ylabel(("Generalized Speeds [Rad/s]"))
plt.title("System Response")

plt.show()
