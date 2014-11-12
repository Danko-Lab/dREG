BWPATH:=$(shell Rscript bigWigpath.R)

CC:=$(shell R CMD config CC)

LDFLAGS:=$(shell R CMD config LDFLAGS)
LIBS:=${shell pkg-config openssl --libs} ${LDFLAGS}

CFLAGS:=-c -g -Wall ${LIBS}
SOURCES= read_genomic_data.c RegisterRRoutines.c
OBJECTS=$(SOURCES:.c=.o)

MFLGS := "CFLAGS=-I${BWPATH}/include ${CFLAGS}"
# define a space
empty:=
space:= $(empty) $(empty)
MFLGSESC := $(subst $(space),\$(space),$(MFLGS))

#.c.o:
#	$(CC) $(CFLAGS) $< -o $@

featureDetector.so:
	MAKEFLAGS=${MFLGSESC} R CMD SHLIB -o $@ ${SOURCES} ${BWPATH}/libs/bigwiglib.a ${LIBS}

clean:
	rm -f *.o *.so
