BWPATH:=$(shell ${R_HOME}/bin/Rscript bigWigpath.R)

CC:=$(shell ${R_HOME}/bin/R CMD config CC)

LDFLAGS:=$(shell ${R_HOME}/bin/R CMD config LDFLAGS)
LIBS:=${shell pkg-config openssl --libs} ${LDFLAGS}

CFLAGS:=-c -g -Wall ${LIBS}
SOURCES= read_genomic_data.c
OBJECTS=$(SOURCES:.c=.o)

MFLGS := "CFLAGS=-I${BWPATH}/include ${CFLAGS}"
# define a space
empty:=
space:= $(empty) $(empty)
MFLGSESC := $(subst $(space),\$(space),$(MFLGS))

#.c.o:
#	$(CC) $(CFLAGS) $< -o $@

dREG.so:
	MAKEFLAGS=${MFLGSESC} ${R_HOME}/bin/R CMD SHLIB -o $@ ${SOURCES} ${BWPATH}/lib/bigwiglib.a ${LIBS}

clean:
	rm -f *.o *.so

