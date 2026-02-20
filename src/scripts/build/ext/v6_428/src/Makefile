#
# simple pythia6 makefile
#
UNAME = $(shell uname)
ifeq "$(UNAME)" "Linux"
    AR=ar
    F77=gfortran
    FFLAG= -O -fno-second-underscore -fPIC 
    CPP = gcc -E
    CPPFLG= -C -P
endif
ifeq "$(UNAME)" "Darwin"
    AR=ar
    F77=gfortran
    FFLAG= -O -fno-second-underscore -fPIC 
    CPP = cc -E
    CPPFLG= -C -P
endif
LIB = ../lib
CPPFLGS = $(CPPFLG) -D$(UNAME)

FOBJSALL = $(patsubst %.f,%.o,$(wildcard *.f))            $(patsubst %.F,%.o,$(wildcard *.F))

# remove the "pdfdum.o" as we don't want that going into the .a library
FOBJS = $(filter-out pdfdum.o,$(FOBJSALL))

#------------------------------------------

all: $(LIB)/liblund.a $(LIB)/pydata.o

$(LIB)/liblund.a: $(FOBJS)
	$(AR) -urs $(LIB)/liblund.a $(FOBJS)

$(LIB)/pydata.o: pydata.o
	cp -p pydata.o $(LIB)/pydata.o

$(LIB)/pdfdum.o: pdfdum.o
	cp -p pdfdum.o $(LIB)/pdfdum.o

clean:
	rm -f *.o

veryclean:
	rm -f $(LIB)/liblund.a $(LIB)/pydata.o $(LIB)/pdfdum.o
	rm -f *.o

#------------------------------------------

.SUFFIXES : .o .f .F

.f.o:
	$(F77) $(FFLAG) -c $<

.F.o:
	$(F77) $(FFLAG) -c $<

