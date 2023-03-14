CC=gfortran

CFLAGS=-c 

LDFLAGS= -L$(LIBRARY_PATH)

all: gipl

gipl: gipl_mods.o gipl.o
	$(CC) gipl_mods.o gipl.o -o gipl

gipl_mods.o: gipl_mods.f90
	$(CC) $(LDFLAGS) $(CFLAGS) gipl_mods.f90

gipl.o: gipl.f90
	$(CC) $(LDFLAGS) $(CFLAGS) gipl.f90

clean:
	rm *o gipl
