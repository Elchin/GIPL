CC=gfortran

CFLAGS=-c 

all: gipl

gipl: gipl_mods.o gipl.o
	$(CC) gipl_mods.o gipl.o -o gipl

gipl_mods.o: gipl_mods.f90
	$(CC) $(CFLAGS) gipl_mods.f90

gipl.o: gipl.f90
	$(CC) $(CFLAGS) gipl.f90

clean:
	rm *o gipl
