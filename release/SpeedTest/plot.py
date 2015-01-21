import numpy as np
import math
import matplotlib.pyplot as plt
class Body:
	pass

with open('x_axis.mtx') as file:
    x_axis = [[float(digit) for digit in line.split()] for line in file]

with open('y_axis.mtx') as file:
    y_axis = [[float(digit) for digit in line.split()] for line in file]

Y=np.array(y_axis,dtype=np.float64)
X=np.array(x_axis,dtype=np.float64)

# Solution Plot
#--------------------
plt.plot(X[0,:],Y[0,:],X[0,:],Y[1,:],X[0,:],Y[2,:],X[0,:],Y[3,:])
plt.xlabel("Number of Bodies")
plt.ylabel("Time [ms]")
plt.title("Speed Test Results")
plt.legend(['Only OpenMP','1 Assembly on GPU', '3 Assemlblies on GPU', '6 Assemblies on GPU'], 'best')
plt.show()

