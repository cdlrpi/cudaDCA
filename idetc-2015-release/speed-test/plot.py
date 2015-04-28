import numpy as np
import math
import matplotlib.pyplot as plt
class Body:
	pass

with open('n.mtx') as file:
    n = np.array([[float(digit) for digit in line.split()] for line in file],dtype=np.float64)
with open('t.mtx') as file:
    t = np.array([[float(digit) for digit in line.split()] for line in file],dtype=np.float64)


# Solution Plot
#--------------------
plt.plot(n[0,:],t[0,:],n[0,:],t[1,:],n[0,:],t[2,:],n[0,:],t[3,:])
plt.xlabel("Number of Bodies")
plt.ylabel("Time [ms]")
plt.title("Speed Test Results")
plt.legend(['Only OpenMP','1 Assembly on GPU', '3 Assemlblies on GPU', '6 Assemblies on GPU'], 'best')
plt.show()

