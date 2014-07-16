.SUFFIXES: .cu .o

# Build tools
NVCC = /usr/local/cuda-5.5/bin/nvcc

# here are all the objects 
OBJS = cudaDCAdriver.o funct_bin/DCAfuncts.o funct_bin/Functs.o funct_bin/Assemble.o funct_bin/solveBCs.o funct_bin/Disassemble.o funct_bin/RK45.o funct_bin/npy.o funct_bin/cudaDCA.o h_code/Assemble_setup.o h_code/Disassemble_setup.o h_code/Init_setup.o funct_bin/Initialize.cu funct_bin/updatezetas.cu

#rule to change .cu into .o
.cu.o:
	$(NVCC) -arch=compute_30 -c $< -o $@

# make and compile
go: $(OBJS)
	$(NVCC) $(OBJS) -o go


#all gpu functions are included instead of compiled as separate objects because nvcc
#returns an error when calls on the gpu are made to a separate file.
cudaDCAdriver.o:  funct_bin/npy.h d_code/deviceFuncts.h d_code/deviceAssemble.h d_code/deviceDisassemble.h d_code/deviceInitialize.h
funct_bin/DCAfuncts.o: 
funct_bin/cudaDCA.o: 
funct_bin/Functs.o: 
funct_bin/Assemble.o: 
funct_bin/solveBCs.o: 
funct_bin/Disassemble.o: 
funct_bin/RK45.o: 
funct_bin/npy.o: funct_bin/npy.h
h_code/Assemble_setup.o: 
h_code/Disassemble_setup.o: 
h_code/Init_setup.o:

#clean command to remove all temporary and unnecessary files
clean:
	rm -f *.o *~ funct_bin/*~ funct_bin/*.o *.mtx *.npy go h_code/*.o h_code/*~  core


