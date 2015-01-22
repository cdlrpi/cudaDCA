Parallelized DCA Simulation README file

The RunSim folder contains a version of DCA used to simulate an n-link pendulum.
/--------------------------------------------------------------------------------------------------------------------------------/

Folder Navigation:
------------------
All files are foind in the SRC folder.  Inside this folder is the main function "Main.cu". 
Inside the SRC folder there are 3 folders: Host_Code, Device_Code, and DCA_Functions.
Host_Code contains the files that govern communication between the GPU and CPU.
Device_Code contains the GPU code files.
DCA_Functions contain all CPU code that governs the DCA process.

Operation:
----------
In order for this program to work properly, CUDA must be properly installed and configured.
If CUDA is installed, follow the instructions below to output matrix data of the test.

1) Navigate to the RunSim folder in the terminal
2) Enter 'make' and wait for the program to be compiled
3) Enter './runSim' and wait for the program to finish running

This will create 2 files named 'Body_Values.mtx' and 'Solution.mtx'

In order to plot these files, python with the matplotlib package must be installed.
If python and matplotlib are properly installed and configured:

1)Enter 'python plot.py' to produce a plot of the results

Output File Form:
-----------------
The output file 'Body_Values.mtx' is of the following form:

	The first number of the first line is the timestep size in seconds.
	All other numbers on the first line are a list of the mass of each body in order starting from the link connected to the inertial frame
	The first number of the second line is the final time in seconds.
	All other numbers on the second line are a list of the length of each link in meters

The output file 'Solution.mtx is of the following form:

	Each line contains a list of the generalized coordinates of each joint, followed by the generalized speeds for the corresponding timestep
	Example (3 links):
		line1: q1 q2 q3 qdot1 qdot2 qdot3
		line2: q1 q2 q3 qdot1 qdot2 qdot3

