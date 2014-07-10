import numpy as np
import math
import matplotlib.pyplot as plt
class Body:
	pass



#Y = np.load('mat.npy')


with open('graph_gpucpu1.mtx') as file:
    time3 = [[float(digit) for digit in line.split()] for line in file]
with open('graph1_init_on_cpu.mtx') as file:
    time4 = [[float(digit) for digit in line.split()] for line in file]
with open('graph2_init_on_cpu.mtx') as file:
    time5 = [[float(digit) for digit in line.split()] for line in file]
with open('graph_gpucpu2.mtx') as file:
    time2 = [[float(digit) for digit in line.split()] for line in file]
with open('numbods.mtx') as file:
    YY = [[float(digit) for digit in line.split()] for line in file]


Time3=np.array(time3,dtype=np.float64)
Time4=np.array(time4,dtype=np.float64)
Time5=np.array(time5,dtype=np.float64)
Time2=np.array(time2,dtype=np.float64)
#print(Vals)
numbods=np.array(YY,dtype=np.float64)

plt.plot(numbods[0,:],Time3[0,:],numbods[0,:],Time3[1,:],numbods[0,:],Time3[2,:],numbods[0,:],Time3[3,:],numbods[0,:],Time2[0,:],numbods[0,:],Time2[1,:],numbods[0,:],Time2[2,:],numbods[0,:],Time2[3,:], numbods[0,:],Time4[0,:])
plt.ylabel("Time [s]")
plt.xlabel("number of bodies")
plt.title("graph")
plt.legend(['all cpu', '1 gpu assembly','2 gpu assembly','3 gpu assembly','4 gpu assembly','5 gpu assembly','6 gpu assembly','7 gpu assembly'],'best')
plt.show()
plt.plot(numbods[0,:],Time4[0,:],numbods[0,:],Time4[1,:],numbods[0,:],Time4[2,:],numbods[0,:],Time4[3,:],numbods[0,:],Time5[0,:],numbods[0,:],Time5[1,:],numbods[0,:],Time5[2,:],numbods[0,:],Time5[3,:], numbods[0,:], Time5[4,:], numbods[0,:], Time5[5,:], numbods[0,:], Time5[6,:], numbods[0,:], Time5[7,:])
plt.ylabel("Time [s]")
plt.xlabel("number of bodies")
plt.title("graph")
plt.legend(['all cpu', '1 gpu assembly','2 gpu assembly','3 gpu assembly','4 gpu assembly','5 gpu assembly','6 gpu assembly','7 gpu assembly','8 gpu assembly','9 gpu assembly','10 gpu assembly','11 gpu assembly'],'best')
plt.show()
