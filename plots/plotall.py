import numpy as np
import math
import matplotlib.pyplot as plt
class Body:
	pass



#Y = np.load('mat.npy')


with open('graph_gpucpu1.mtx') as file:
    time3 = [[float(digit) for digit in line.split()] for line in file]

with open('numbods.mtx') as file:
    YY = [[float(digit) for digit in line.split()] for line in file]


Time3=np.array(time3,dtype=np.float64)
#print(Vals)
numbods=np.array(YY,dtype=np.float64)

plt.plot(numbods[0,:],Time3[0,:],numbods[0,:],Time3[1,:],numbods[0,:],Time3[2,:],numbods[0,:],Time3[3,:])
plt.xlabel("Time [s]")
plt.ylabel("number of bodies")
plt.title("graph")
plt.legend(['all cpu', '1 gpu assembly','2 gpu assembly','3 gpu assembly'],'best')
plt.show()
