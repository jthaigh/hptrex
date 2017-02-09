#=========================================================
# Makefile
#=========================================================

#CC	= g++
COPTS	= -Wall -g -fPIC -DLINUX -O2 -std=c++11 -Wl,-rpath,$(shell brew --prefix)/lib

#LIBS	= -L/usr/local/lib -lftd2xx

OBJ = TTPCAStar.o TTPCHitPad.o TTPCLayout.o TTPCOrderedVolGroup.o TTPCPathVolume.o TTPCTRExPatAlgorithm.o TTPCTRExPatSubAlgorithm.o TTPCUnitVolume.o TTPCVolGroup.o TTPCVolGroupMan.o GasTPCDataLib.o TSimLoader.o GasTPCDataLibDict.o TTRExPattern.o TCCDSimLoader.o TTRExPatternDict.o

ROOT_LIBS   = $(shell root-config --libs --glibs)
ROOT_CFLAGS = $(shell root-config --cflags)

#=========================================================
# Compile
#=========================================================
all: RunTRExRecon_CCD.exe RunTRExRecon.exe
 
RunTRExRecon_CCD.exe: libpatrec.so RunTRExRecon_CCD.cxx
	$(CC) $(COPTS) $(ROOT_CFLAGS) $(ROOT_LIBS) $(INCLUDE) -I. -L. -lpatrec -o RunTRExRecon_CCD.exe RunTRExRecon_CCD.cxx

RunTRExRecon.exe: libpatrec.so RunTRExRecon.cxx
	$(CC) $(COPTS) $(ROOT_CFLAGS) $(ROOT_LIBS) $(INCLUDE) -I. -L. -lpatrec -o RunTRExRecon.exe RunTRExRecon.cxx

libpatrec.so: $(OBJ)
	$(CC) $(COPTS) -shared -o libpatrec.so $(OBJ)

#obj: $(OBJ)

%.o :	%.cxx
	$(CC) -c $(COPTS) $(ROOT_CFLAGS) $(INCLUDE) -I. -o $@ $<

GasTPCDataLibDict.cxx:
	rootcint GasTPCDataLibDict.cxx GasTPCDataLib.hxx LinkDef.hh

#TTRExPattern.o: TTRExPattern.hxx

TTRExPatternDict.cxx: 
	rootcint TTRExPatternDict.cxx  TTRExPattern.hxx OutLink.hh

clean:
	rm *.o *.so *.exe