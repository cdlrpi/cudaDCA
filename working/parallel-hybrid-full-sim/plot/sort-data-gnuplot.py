import numpy as np

with open('../parallel/parallelBodies.mtx') as file:
    n = np.array([[float(digit) for digit in line.split()] for line in file],dtype=np.float64)
with open('../parallel/parallelTimes.mtx') as file:
    tp = np.array([[float(digit) for digit in line.split()] for line in file],dtype=np.float64)
with open('../serial/serialTimes.mtx') as file:
    ts = np.array([[float(digit) for digit in line.split()] for line in file],dtype=np.float64)


# Output data in GNU-Plot friendly format
#--------------------
data = np.hstack((np.reshape(n[0],(n.shape[1],1)),ts.T,tp.T))
np.savetxt('data.mat', data, fmt='%3.16f')
