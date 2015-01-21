.SUFFIXES: .cu .o

# Build tools
NVCC = /usr/local/cuda-5.5/bin/nvcc

# here are all the objects 
OBJS = Assemble_setup.o Disassemble_setup.o Assemble.o Disassemble.o findCutoff.o DCAfuncts.o 

#rule to change .cu into .o
.cu.o:
	$(NVCC) -arch=compute_30 -c $< -o $@

# make and compile
go: $(OBJS)
	$(NVCC) $(OBJS) -o go


#all gpu functions are included instead of compiled as separate objects because nvcc
#returns an error when calls on the gpu are made to a separate file.
findCutoff.o: deviceAssemble.h deviceDisassemble.h deviceFuncts.h
Assemble_setup.o:
Disassemble_setup.o:
Assemble.o:
Disassemble.o:

#clean command to remove all temporary and unnecessary files
clean:
	rm -f *.o *~ funct_bin/*~ funct_bin/*.o *.mtx *.npy go h_code/*.o h_code/*~  core


