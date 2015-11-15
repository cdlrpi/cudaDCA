import numpy as np
import math
import matplotlib.pyplot as plt

with open('../parallel-hybrid-speed-test/n.mtx') as file:
    n = np.array([[float(digit) for digit in line.split()] for line in file],dtype=np.float64)
with open('../parallel-hybrid-speed-test/t.mtx') as file:
    tp = np.array([[float(digit) for digit in line.split()] for line in file],dtype=np.float64)                   
with open('../serial-speed-test/t.mtx') as file:
    ts = np.array([[float(digit) for digit in line.split()] for line in file],dtype=np.float64)          


# Solution Plot
#--------------------
plt.plot(n[0,:],ts[0,:],n[0,:],tp[0,:],n[0,:],tp[1,:],n[0,:],tp[2,:],n[0,:],tp[3,:])
plt.xlabel("Number of Bodies")
plt.ylabel("Time [ms]")
plt.title("Speed Test Results")
plt.legend(['Serial','Only OpenMP','1 Assembly on GPU', '3 Assemlblies on GPU', '6 Assemblies on GPU'], loc='best')
plt.show()

