CC     = mpicc
CFLAGS = -O3
EFLAGS = 
EFILE  = PyroNoiseA
LIBS   =  -lm
OBJS   = PyroNoiseA.o 

$(EFILE) : $(OBJS)
	@echo "linking..."
	$(CC) $(EFLAGS) -o $(EFILE) $(OBJS) $(LIBS)

$(OBJS) : PyroNoiseA.c PyroNoiseA.h
	$(CC) $(CFLAGS) -c $*.c 

clean:
	rm -rf *.o PyroNoiseA
