OBJS= convolutionMPI.o 
BIN= convolutionMPI


CC=mpicc
CXX=mpiCC
NVCC=nvcc
NVFLAGS= -ccbin=$(CXX)
NVLDFLAGS=


all: $(OBJS)
	$(NVCC) $(NVFLAGS) $(NVLDFLAGS) $(OBJS) -o $(BIN)

%.o: %.c
	$(CC) $(CFLAGS)  -c $< -o $(notdir $@ ) 

%.o: %.cpp
	$(CXX) $(CFLAGS)  -c $< -o $(notdir $@ ) 

%.o: %.cu
	$(NVCC) $(NVFLAGS)  -c $< -o $(notdir $@ )

clean:
	rm $(OBJS) $(BIN)
