CC=nvcc
CFLAGS=-Wall -O3 
% : %.cc	
	    $(CC) $(CFLAGS) $< -o $@
