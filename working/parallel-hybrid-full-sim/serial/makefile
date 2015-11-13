.SUFFIXES: .cu .o

# Build tools
NVCC = /usr/local/cuda/bin/nvcc

# here are all the objects 
OBJS = SRC/Main.o SRC/DCA_Functions/DCAfuncts.o SRC/DCA_Functions/Functs.o SRC/DCA_Functions/Assemble.o SRC/DCA_Functions/solveBCs.o SRC/DCA_Functions/Disassemble.o SRC/DCA_Functions/RK45.o  SRC/DCA_Functions/hybridDCA.o SRC/Host_Code/Assemble_setup.o SRC/Host_Code/Disassemble_setup.o SRC/Host_Code/Init_setup.o SRC/DCA_Functions/Initialize.cu SRC/DCA_Functions/updatezetas.cu SRC/DCA_Functions/findCutoff.cu

#rule to change .cu into .o
.cu.o:
	$(NVCC)  -c $< -o $@

.cpp.o:
	$(NVCC) -Xcompiler -fopenmp -c $< -o $@

# make and compile
go: $(OBJS)
	$(NVCC) $(OBJS) -lgomp -o  runSim


#all gpu functions are included instead of compiled as separate objects because nvcc
#returns an error when calls on the gpu are made to a separate file.
Main.o: SRC/Device_Code/deviceFuncts.h SRC/Device_Code/deviceAssemble.h SRC/Device_Code/deviceDisassemble.h SRC/Device_Code/deviceInitialize.h

#clean command to remove all temporary and unnecessary files
clean:
	rm -f *.o *~ SRC/*~ SRC/*.o SRC/DCA_Functions/*~ SRC/DCA_Functions/*.o *.mtx runSim SRC/Host_Code/*.o SRC/Host_Code/*~  core


