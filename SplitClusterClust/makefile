CC     = gcc
CFLAGS = -O3
EFLAGS =  
EFILE  = SplitClusterClust
LIBS   = -lm 
OBJS   = SplitClusterClust.o

$(EFILE) : $(OBJS)
	@echo "linking..."
	$(CC) $(EFLAGS) -o $(EFILE) $(OBJS) $(LIBS)

$(OBJS) : SplitClusterClust.c SplitClusterClust.h
	$(CC) $(CFLAGS) -c $*.c 

clean:
	rm -rf *.o SplitClusterClust
