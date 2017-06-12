#=========================================================
# Makefile
#=========================================================

#CC	= g++
COPTS	= -Wall -g -fPIC -DLINUX -O2 -std=c++11 -Wl,-rpath,$(shell brew --prefix)/lib

#LIBS	= -L/usr/local/lib -lftd2xx

OBJ = TTPCAStar.o TTPCHitPad.o TTPCLayout.o TTPCOrderedVolGroup.o TTPCPathVolume.o TTPCTRExPatAlgorithm.o TTPCTRExPatSubAlgorithm.o TTPCUnitVolume.o TTPCVolGroup.o TTPCVolGroupMan.o GasTPCDataLib.o TSimLoaderCCD.o TSimLoaderND.o GasTPCDataLibDict.o TTRExPattern.o TEventDisplay.o ioDict.o TTRExHVCluster.o TTPCHelixPropagator.o TTPCLikelihoodMatch.o TTPCLikelihoodMerge.o TTPCLinearMerge.o TTPCLikFitPath.o TTPCQLikelihood.o TTPCSeeding.o TTPCTracking.o TTPCUtils.o TTrueHit.o TTRExPath.o TTRExJunction.o TTRExPIDAlgorithm.o TTrueTrack.o

ROOT_LIBS   = $(shell root-config --libs --glibs)
ROOT_CFLAGS = $(shell root-config --cflags)

#=========================================================
# Compile
#=========================================================

all: libevtrecord.so RunTRExReconCCD.exe RunTRExReconND.exe

libevtrecord.so: evtrec
	make -C evtrec
	cp evtrec/libevtrecord.so ./

#RunTRExRecon.exe: libpatrec.so libevtrecord.so RunTRExRecon.cxx
#	$(CC) $(COPTS) $(ROOT_CFLAGS) $(ROOT_LIBS) $(INCLUDE) -I. -L. -levtrecord -lpatrec -lMinuit -o RunTRExRecon.exe RunTRExRecon.cxx 

RunTRExReconCCD.exe: libpatrec.so RunTRExReconCCD.cxx
	$(CC) $(COPTS) $(ROOT_CFLAGS) $(ROOT_LIBS) $(INCLUDE) -I. -L. -levtrecord -lpatrec -lMinuit -o RunTRExReconCCD.exe RunTRExReconCCD.cxx

RunTRExReconND.exe: libpatrec.so libevtrecord.so RunTRExReconND.cxx 
	$(CC) $(COPTS) $(ROOT_CFLAGS) $(ROOT_LIBS) $(INCLUDE) -I. -L. -levtrecord -lpatrec -lMinuit -o RunTRExReconND.exe RunTRExReconND.cxx

libpatrec.so: $(OBJ)
	$(CC) $(COPTS) -shared -o libpatrec.so $(OBJ)

#obj: $(OBJ)

%.o :	%.cxx
	$(CC) -c $(COPTS) $(ROOT_CFLAGS) $(INCLUDE) -I. -o $@ $<

GasTPCDataLibDict.cxx:
	rootcint GasTPCDataLibDict.cxx GasTPCDataLib.hxx LinkDef.hh

ioDict.cxx: 
	rootcint ioDict.cxx  TTRExPattern.hxx TTPCHitPad.hxx ioLinkDef.hh

clean:
	make -C evtrec clean	
	rm -f *.pcm *.o *.so *.exe GasTPCDataLibDict.* ioDict.* ioDict_rdict.pcm
	