.SUFFIXES: .cu .o

# Build tools
NVCC = /usr/local/cuda/bin/nvcc

# here are all the objects 
OBJS = src/main.o src/dca-functs/dca-functs.o src/dca-functs/functs.o src/dca-functs/assemble.o src/dca-functs/solveBCs.o src/dca-functs/disassemble.o src/dca-functs/rk45.o  src/dca-functs/cudaDCA.o src/h-code/assemble-setup.o src/h-code/disassemble-setup.o src/h-code/Init_setup.o src/dca-functs/initialize.cu src/dca-functs/updatezetas.cu

#rule to change .cu into .o
.cu.o:
	$(NVCC)  -c $< -o $@

.cpp.o:
	$(NVCC) -Xcompiler -fopenmp -c $< -o $@

# make and compile
go: $(OBJS)
	$(NVCC) $(OBJS) -lgomp -o  runSpeedTest


#all gpu functions are included instead of compiled as separate objects because nvcc
#returns an error when calls on the gpu are made to a separate file.
main.o: src/d-code/deviceFuncts.h src/d-code/deviceAssemble.h src/d-code/deviceDisassemble.h src/d-code/deviceInitialize.h

#clean command to remove all temporary and unnecessary files
clean:
	rm -f *.o *~ src/*~ src/*.o src/dca-functs/*~ src/dca-functs/*.o *.mtx runSpeedTest src/h-code/*.o src/h-code/*~  core


