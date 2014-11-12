KENTINC:=kent/src/inc/
KENTOBJS:=kent/src/lib/*.o

all: bigwiglib.a

bigwiglib.a: bigwiglib.h bigwiglib.c
	gcc -Wall -g -fPIC -c bigwiglib.c -o bigwiglib.o -I$(KENTINC)
	ar rcs $@ bigwiglib.o $(KENTOBJS)

clean:
	rm -f *.o *.a
