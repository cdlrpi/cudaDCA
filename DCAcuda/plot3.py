

import numpy as np
import math
import matplotlib.pyplot as plt
class Body:
	pass



#Y = np.load('mat.npy')
with open('graph_gpu.mtx') as file:
    time = [[float(digit) for digit in line.split()] for line in file]

with open('numbods.mtx') as file:
    YY = [[float(digit) for digit in line.split()] for line in file]
Time=np.array(time,dtype=np.float64)
#print(Vals)
numbods=np.array(YY,dtype=np.float64)
print(Time)
print(numbods)
plt.plot(Time[0,:],numbods[0,:])
plt.xlabel("Time [s]")
plt.ylabel("number of bodies")
plt.title("graph")
plt.show()
