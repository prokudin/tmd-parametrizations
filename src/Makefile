#
# Author: Alexei Prokudin <prokudin@jlab.org>
#
CPP = @g++
FF  = @gfortran

FFLAGS = 

ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)
ROOTGLIBS  := $(shell root-config --glibs)

LDFLAGS = -O
CFLAGS += $(ROOTCFLAGS)

# In order for gfortran to work I need to include this
FLIBS =  -L/usr/local/gfortran/lib/ -lgfortran

LIB = $(ROOTLIBS) $(FLIBS)

GRV98    = ./grv98.o
TRANS_LO = ./transv_pdf.o
FRAGMENT = fDSS.o pkhff.o kkp.o polin2.o polint.o grille_had_charged.o locate.o dlib.o akk.o

OPT = -c -O3 -Wall -I./

OBJECTS = parameters.o unpolarised.o sivers.o transversity.o collins.o stfunctions.o \
	$(GRV98) $(TRANS_LO) $(FRAGMENT)


#We create here our library for particles
library: $(OBJECTS)
	ar rvs libTMDS.a $(OBJECTS) 

parameters.o:	parameters.cpp 
	$(CPP) $(OPT) -o $@ parameters.cpp
	@echo "..................done Parameters."

unpolarised.o:	unpolarised.cpp 
	$(CPP) $(OPT) -o $@ unpolarised.cpp
	@echo "..................done Unpolarised."

sivers.o:	sivers.cpp 
	$(CPP) $(OPT) $(CFLAGS) -o $@ sivers.cpp
	@echo "..................done Sivers."

transversity.o:	transversity.cpp 
	$(CPP) $(OPT) $(CFLAGS) -o $@ transversity.cpp
	@echo "..................done Transversity."

collins.o:	collins.cpp 
	$(CPP) $(OPT) $(CFLAGS) -o $@ collins.cpp
	@echo "..................done Collins."

stfunctions.o:	stfunctions.cpp 
	$(CPP) $(OPT) $(CFLAGS) -o $@ stfunctions.cpp
	@echo "..................done Structure Functions."

grv98.o: grv98.f
	$(FF) -c grv98.f

toy.o: toy.f
	$(FF) -c toy.f

transv_pdf.o: transv_pdf.f
	$(FF) -c transv_pdf.f

fDSS.o:  fDSS.f
	$(FF) -c fDSS.f

dlib.o:  dlib.f
	$(FF) -c dlib.f

pkhff.o: pkhff.f 
	$(FF) -c pkhff.f

kkp.o: kkp.f
	$(FF) -c kkp.f

polin2.o: polin2.f 
	$(FF) -c polin2.f

polint.o: polint.f 
	$(FF) -c polint.f

grille_had_charged.o:  grille_had_charged.f
	$(FF) -c grille_had_charged.f

locate.o: locate.f 
	$(FF) -c locate.f
 
akk.o: akk.f
	$(FF) -c akk.f

toy.run: toy.o $(FRAGMENT)
	$(FF) $(FFLAGS) -o $@ $+

make all:
	make toy.run

clean:
	rm *.o
	rm *.run
	rm *.a