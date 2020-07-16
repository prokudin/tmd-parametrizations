#
# Author: Alexei Prokudin <prokudin@jlab.org>
#
CPP = @g++
FF  = @gfortran


ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)
ROOTGLIBS  := $(shell root-config --glibs)

LDFLAGS = -O
CFLAGS += $(ROOTCFLAGS)

# In order for gfortran to work I need to include this
FLIBS =  -L/usr/local/gfortran/lib/ -lgfortran

#tmds
LIB_TMDS = -L./src -lTMDS

LIBS = $(ROOTLIBS) $(LIB_TMDS) $(FLIBS)

OPT = -c -O3 -Wall -I./

test.exe:	test.o
	$(CPP) -o test.exe test.o $(LIBS)

test.o:	test.cpp 
	$(CPP) $(OPT) $(CFLAGS) -o $@ test.cpp
	@echo "..................done Test."

stfunctionstest.exe:	stfunctionstest.o
	$(CPP) -o stfunctiontest.exe stfunctionstest.o $(LIBS)

stfunctionstest.o:	stfunctionstest.cpp 
	$(CPP) $(OPT) $(CFLAGS) -o $@ stfunctionstest.cpp
	@echo "..................done Structure Functions test."

make all:
	$(MAKE) --directory=src
	make test.exe
	make stfunctionstest.exe

clean:
	rm *.o
	rm *.exe
	$(MAKE) --directory=src clean
