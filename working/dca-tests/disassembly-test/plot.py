import numpy as np
import math
import matplotlib.pyplot as plt
class Body:
	pass



#Y = np.load('mat.npy')



with open('cpudisassembletime2.mtx') as file:
    time4 = [[float(digit) for digit in line.split()] for line in file]
with open('gpudisassembletime2.mtx') as file:
    time5 = [[float(digit) for digit in line.split()] for line in file]
with open('numbods.mtx') as file:
    YY = [[float(digit) for digit in line.split()] for line in file]

with open('gpudisassembletime.mtx') as file:
    time6 = [[float(digit) for digit in line.split()] for line in file]

Time4=np.array(time4,dtype=np.float64)
Time5=np.array(time5,dtype=np.float64)
Time6=np.array(time6,dtype=np.float64)

#print(Vals)
numbods=np.array(YY,dtype=np.float64)

#plt.plot(numbods[0,:],Time4[0,:],numbods[0,:],Time4[1,:],numbods[0,:],Time4[2,:],numbods[0,:],Time4[3,:],numbods[0,:],Time4[4,:],numbods[0,:],Time5[0,:],numbods[0,:],Time5[1,:],numbods[0,:],Time5[2,:],numbods[0,:],Time5[3,:], numbods[0,:], Time5[4,:], numbods[0,:], Time5[5,:], numbods[0,:], Time5[6,:])
plt.plot(numbods[0,:],Time4[0,:],numbods[0,:],Time5[0,:], numbods[0,:],Time6[0,:])
plt.ylabel("Time [ms]")
plt.xlabel("number of bodies")
plt.title("gpu vs cpu disassembly")
plt.legend(['cpu disassembly', 'gpu disassembly with transfer','gpu disassembly no transfer'],'best')
plt.show()
