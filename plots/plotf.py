import numpy as np
import math
import matplotlib.pyplot as plt
class Body:
	pass



#Y = np.load('mat.npy')



with open('4kcudadcalab2.mtx') as file:
    cuda = [[float(digit) for digit in line.split()] for line in file]
with open('numbods4klab2.mtx') as file:
    YY = [[float(digit) for digit in line.split()] for line in file]
with open('ompdcalab2.mtx') as file:
    omp = [[float(digit) for digit in line.split()] for line in file]


CUDA=np.array(cuda,dtype=np.float64)


#print(Vals)
numbods=np.array(YY,dtype=np.float64)

OMP = np.array(omp,dtype=np.float64)
#plt.plot(numbods[0,:],Time4[0,:],numbods[0,:],Time4[1,:],numbods[0,:],Time4[2,:],numbods[0,:],Time4[3,:],numbods[0,:],Time4[4,:],numbods[0,:],Time5[0,:],numbods[0,:],Time5[1,:],numbods[0,:],Time5[2,:],numbods[0,:],Time5[3,:], numbods[0,:], Time5[4,:], numbods[0,:], Time5[5,:], numbods[0,:], Time5[6,:])
plt.plot(numbods[0,:],CUDA[0,:],numbods[0,:],CUDA[1,:],numbods[0,:],CUDA[2,:],numbods[0,:],CUDA[3,:],numbods[0,:],OMP[0,:],numbods[0,:],OMP[1,:],numbods[0,:],OMP[2,:],numbods[0,:],OMP[3,:])
#plt.plot(numbods[0,:],OMP[0,:],numbods[0,:],OMP[1,:],numbods[0,:],OMP[2,:],numbods[0,:],OMP[3,:])
plt.ylabel("Time [ms]")
plt.xlabel("number of bodies")
plt.title("Time to Solve Using DCA: CUDA and openMP Comparison using lab machine")
plt.legend(['Single CPU', '1 lvl on gpu','3 lvl on gpu','6 lvl on gpu','8 Core openMP','1 lvl on gpu w/ omp','3 lvl on gpu w/ omp','6 lvl on gpu w/ omp'],'best')
plt.show()


