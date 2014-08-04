import numpy as np
import math
import matplotlib.pyplot as plt
class Body:
	pass



#Y = np.load('mat.npy')



with open('graph_hdca1-1f.mtx') as file:
    cuda = [[float(digit) for digit in line.split()] for line in file]
with open('numbods1-1f.mtx') as file:
    YY = [[float(digit) for digit in line.split()] for line in file]
with open('graphomphyb.mtx') as file:
    omp = [[float(digit) for digit in line.split()] for line in file]

print(time4)
cuda=np.array(CUDA,dtype=np.float64)


#print(Vals)
numbods=np.array(YY,dtype=np.float64)

omp = np.array(OMP,dtype=np.float64)
#plt.plot(numbods[0,:],Time4[0,:],numbods[0,:],Time4[1,:],numbods[0,:],Time4[2,:],numbods[0,:],Time4[3,:],numbods[0,:],Time4[4,:],numbods[0,:],Time5[0,:],numbods[0,:],Time5[1,:],numbods[0,:],Time5[2,:],numbods[0,:],Time5[3,:], numbods[0,:], Time5[4,:], numbods[0,:], Time5[5,:], numbods[0,:], Time5[6,:])
plt.plot(numbods[0,:],CUDA[0,:],numbods[0,:],CUDA[1,:],numbods[0,:],CUDA[3,:],numbods[0,:],CUDA[6,:],numbods[0,:],OMP[0,:],numbods[0,:],OMP[1,:],numbods[0,:],OMP[3,:],numbods[0,:],OMP[6,:])
plt.ylabel("Time [s]")
plt.xlabel("number of bodies")
plt.title("graph")
plt.legend(['all cpu', '1 gpu assembly','2 gpu assembly','3 gpu assembly','4 gpu assembly','5 gpu assembly','6 gpu assembly','7 gpu assembly','8 gpu assembly','9 gpu assembly','10 gpu assembly','11 gpu assembly'],'best')
plt.show()


