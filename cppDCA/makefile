.SUFFIXES: .cu .o

# Build tools
NVCC = /usr/local/cuda-5.5/bin/nvcc

# here are all the objects 
OBJS = DCAdriver.o bin/Initialize.o bin/DCAfuncts.o bin/DCA.o bin/Functs.o bin/Assemble.o bin/solveBCs.o bin/Disassemble.o bin/RK45.o bin/npy.o

#rule to change .cu into .o
.cu.o:
	$(NVCC) -c $< -o $@

# make and compile
go: $(OBJS)
	$(NVCC) $(OBJS) -o go


#all gpu functions are included instead of compiled as separate objects because nvcc
#returns an error when calls on the gpu are made to a separate file.
DCAdriver.o: bin/classes.h bin/npy.h
Initialize.o: classes.h
DCAfuncts.o: classes.h
DCA.o: classes.h
Functs.o: classes.h
Assemble.o: classes.h
solveBCs.o: classes.h
Disassemble.o: classes.h
RK45.o: classes.h
npy.o: npy.h
#RK45.o: classes.h
#clean command to remove all temporary and unnecessary files
clean:
	rm -f *.o *~ bin/*~ bin/*.o *.mtx *.npy go  core


