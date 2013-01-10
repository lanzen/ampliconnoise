CC     = cc
CFLAGS = -O3 
EFLAGS =  
EFILE  = PerseusD
LIBS   = -lm -lgsl -lgslcblas
OBJS   = PerseusD.o 

$(EFILE) : $(OBJS)
	@echo "linking..."
	$(CC) $(EFLAGS) -o $(EFILE) $(OBJS) $(LIBS)

$(OBJS) : PerseusD.c PerseusD.h
	$(CC) $(CFLAGS) -c $*.c 

clean:
	rm -rf *.o PerseusD
