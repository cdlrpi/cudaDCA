

import numpy as np
import math
import matplotlib.pyplot as plt
class Body:
	pass



#Y = np.load('mat.npy')
with open('graph_gpu.mtx') as file:
    time = [[float(digit) for digit in line.split()] for line in file]

with open('graph_cpu.mtx') as file:
    time2 = [[float(digit) for digit in line.split()] for line in file]


with open('numbods.mtx') as file:
    YY = [[float(digit) for digit in line.split()] for line in file]
Time1=np.array(time,dtype=np.float64)
Time2=np.array(time2,dtype=np.float64)
#print(Vals)
numbods=np.array(YY,dtype=np.float64)

plt.plot(Time1[0,:],numbods[0,:],Time2[0,:],numbods[0,:])
plt.xlabel("Time [s]")
plt.ylabel("number of bodies")
plt.title("graph")
plt.legend(['gpu','cpu'],'best')
plt.show()
