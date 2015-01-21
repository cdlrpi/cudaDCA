import numpy as np
import math
import matplotlib.pyplot as plt
class Body:
	pass



#Y = np.load('mat.npy')



with open('graphomphyb.mtx') as file:
    time4 = [[float(digit) for digit in line.split()] for line in file]
with open('numbodsomphyb.mtx') as file:
    YY = [[float(digit) for digit in line.split()] for line in file]


print(time4)
Time4=np.array(time4,dtype=np.float64)


#print(Vals)
numbods=np.array(YY,dtype=np.float64)

#plt.plot(numbods[0,:],Time4[0,:],numbods[0,:],Time4[1,:],numbods[0,:],Time4[2,:],numbods[0,:],Time4[3,:],numbods[0,:],Time4[4,:],numbods[0,:],Time5[0,:],numbods[0,:],Time5[1,:],numbods[0,:],Time5[2,:],numbods[0,:],Time5[3,:], numbods[0,:], Time5[4,:], numbods[0,:], Time5[5,:], numbods[0,:], Time5[6,:])
plt.plot(numbods[0,:],Time4[0,:],numbods[0,:],Time4[1,:],numbods[0,:],Time4[2,:],numbods[0,:],Time4[3,:],numbods[0,:],Time4[4,:],numbods[0,:],Time4[5,:],numbods[0,:],Time4[6,:])
plt.ylabel("Time [s]")
plt.xlabel("number of bodies")
plt.title("graph")
plt.legend(['all cpu', '1 gpu assembly','2 gpu assembly','3 gpu assembly','4 gpu assembly','5 gpu assembly','6 gpu assembly','7 gpu assembly','8 gpu assembly','9 gpu assembly','10 gpu assembly','11 gpu assembly'],'best')
plt.show()


