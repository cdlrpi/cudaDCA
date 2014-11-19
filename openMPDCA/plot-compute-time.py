import numpy as np
import math
import matplotlib.pyplot as plt
class Body:
	pass



#Y = np.load('mat.npy')



with open('time.mtx') as file:
    time = [[float(digit) for digit in line.split()] for line in file]
with open('n-bodies.mtx') as file:
    bodies = [[float(digit) for digit in line.split()] for line in file]

# Convert to numpy arrays of proper size
time=np.array(time,dtype=np.float64)
bodies=np.array(bodies,dtype=np.float64)

plt.plot(bodies[0,:],time[0,:],bodies[0,:],time[1,:],bodies[1,:],time[2,:],bodies[1,:],time[3,:])
plt.ylabel("Time [ms]")
plt.xlabel("number of bodies")
plt.title("Compute Time of DCA using Various Parallel Hybrids")
plt.legend(['All OpenMP', '1 lvl on GPU','6 lvls on GPU','All GPU'],'best')
plt.show()


