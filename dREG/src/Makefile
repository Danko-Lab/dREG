ifeq ($(BIGWIGLIBHOME),)
       BIGWIGLIBHOME:=lib/
endif

ifeq (${MACHTYPE},)
        MACHTYPE=x86_64
endif

CC=gcc
LIBS:=-fPIC -lm -lz -pthread ${shell pkg-config openssl --libs}
CFLAGS:=-fPIC -c -g -Wall ${LIBS}
SOURCES= read_genomic_data.c RegisterRRoutines.c
OBJECTS=$(SOURCES:.c=.o)

MFLGS := "CFLAGS=-I${BIGWIGLIBHOME} -I${BIGWIGLIBHOME}/kent/src/inc ${CFLAGS}"
# define a space
empty:=
space:= $(empty) $(empty)
MFLGSESC := $(subst $(space),\$(space),$(MFLGS))

#.c.o:
#	$(CC) $(CFLAGS) $< -o $@

dREG.so:  
	MAKEFLAGS=${MFLGSESC} R CMD SHLIB -fPIC -o $@ ${SOURCES} ${BIGWIGLIBHOME}/bigwiglib.a ${LIBS}

clean:
	rm -f *.o *.so

