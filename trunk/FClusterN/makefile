CC     = cc
CFLAGS = -O3
EFLAGS =  
EFILE  = FClusterN
LIBS   = -lm 
OBJS   = FClusterN.o 

$(EFILE) : $(OBJS)
	@echo "linking..."
	$(CC) $(EFLAGS) -o $(EFILE) $(OBJS) $(LIBS)

$(OBJS) : FClusterN.c FClusterN.h
	$(CC) $(CFLAGS) -c $*.c 

clean: 
	rm -rf *o FClusterN
