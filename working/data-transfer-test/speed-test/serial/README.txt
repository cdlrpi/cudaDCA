Parallelized DCA Speet Test README file

The SpeedTest folder contains a stripped down version of DCA used strictly to compare the speed of various levels of GPU usage.
/--------------------------------------------------------------------------------------------------------------------------------/

Folder Navigation:
------------------
All files are foind in the SRC folder.  Inside this folder is the main function "Main.cu". 
Inside the SRC folder there are 3 folders: Host_Code, Device_Code, and DCA_Functions.
Host_Code contains the files that govern communication between the GPU and CPU.
Device_Code contains the GPU code files.
DCA_Functions contain all CPU code that governs the DCA process.

Settings:
---------
The current settings run 4 speed tests for a general comparison.
The tests performed are:
	1) All assemblies are done parallelized with openMP on the CPU
	2) 1 assembly is done parallelized on the GPU, all other assemblies are done parallelized with openMP on the CPU
	3) 3 assemblies are done parallelized on the GPU, all other assemblies are done parallelized with openMP on the CPU
	4) 6 assemblies are done parallelized on the GPU, all other assemblies are done parallelized with openMP on the CPU

To alter these settings, open Main.cu and read the comments describing the meaning of the variables and loops.

Operation:
----------
In order for this program to work properly, CUDA must be properly installed and configured.
If CUDA is installed, follow the instructions below to output matrix data of the test.

1) Navigate to the SpeedTest folder in the terminal
2) Enter 'make' and wait for the program to be compiled
3) Enter './runSpeedTest' and wait for the program to finish running

This will create 2 files named 'x_axis.mtx' and 'y_axis.mtx'

In order to plot these files, python with the matplotlib package must be installed.
If python and matplotlib are properly installed and configured:

1)Enter 'python plot.py' to produce a plot of the results

Output File Form:
-----------------
The output file 'x_axis.mtx' is of the following form:

	Each line contains a list of the number of bodies for each speed test data point separated by spaces.
	Example:
		line1: 5  10  15  20  25  30  35  40 ...
		line2: 5  10  15  20  25  30  35  40 ...
	Where line1 corresponds to the first speed test (All assemblies perfomed on the CPU) and each number corresponds to the number of bodies for that data point

The output file 'y_axis.mtx is of the following form:

	Each line contains a list of the time (in ms) to complete a single DCA run for the corresponding number of bodies, separated by spaces.
	Example:
		line1: 0.155802  0.20528  0.249645  0.51929  0.404832  0.414374 ...
		line2: 2.78548  2.07186  2.07452  1.75265  2.04756  1.62085 ...
	Where line1 corresponds to the first speed test (All assemblies performed on the CPU) and each number corresponds to the time of completiong for that data point
