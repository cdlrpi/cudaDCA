.SUFFIXES: .cu .o

# Build tools
NVCC = /usr/local/cuda-5.5/bin/nvcc

# here are all the objects 
OBJS = cudaDCAdriver.o h_code/Disassemble_setup.o h_code/Init_setup.o h_code/Assemble_setup.o h_code/solveBCs.o h_code/cudaDCA.o h_code/RK45.o h_code/npy.o	h_code/Functs.o

#rule to change .cu into .o
.cu.o:
	$(NVCC) -arch=compute_30 -c $< -o $@

# make and compile
go: $(OBJS)
	$(NVCC) $(OBJS) -o go


#all gpu functions are included instead of compiled as separate objects because nvcc
#returns an error when calls on the gpu are made to a separate file.
cudaDCAdriver.o: h_code/classes.h h_code/deviceHeader.h d_code/deviceFuncts.h d_code/deviceAssemble.h d_code/deviceDisassemble.h d_code/deviceInitialize.h h_code/npy.h
h_code/Assemble_setup.o: h_code/classes.h
h_code/Disassemble_setup.o: h_code/classes.h
h_code/Init_setup.o: h_code/classes.h
h_code/solveBCs.o: h_code/classes.h
h_code/cudaDCA.o: h_code/classes.h
h_code/RK45.o: h_code/classes.h
h_code/npy.o: h_code/npy.h
h_code/Functs.o: h_code/classes.h
#RK45.o: classes.h
#clean command to remove all temporary and unnecessary files
clean:
	rm -f *.o *~ h_code/*~ h_code/*.o d_code/*~ d_code/*.o go *.txt core


