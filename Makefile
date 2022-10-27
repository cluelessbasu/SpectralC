CC=g++
CFLAGS=-O3 
% : %.cpp	
	    $(CC) $(CFLAGS) $< -o $@
