import numpy as np

with open('n.mtx') as file:
    n = np.array([[float(digit) for digit in line.split()] for line in file],dtype=np.float64)
with open('t.mtx') as file:
    t = np.array([[float(digit) for digit in line.split()] for line in file],dtype=np.float64)


# Output data in GNU-Plot friendly format
#--------------------
data = np.hstack((np.reshape(n[0],(n.shape[1],1)),t.T))
np.savetxt('data.mat', data, fmt='%3.16f')
