#=========================================================
# Makefile
#=========================================================

#CC	= g++
COPTS	= -Wall -g -fPIC -DLINUX -O2 -std=c++11 -Wl,-rpath,$(shell brew --prefix)/lib

#LIBS	= -L/usr/local/lib -lftd2xx

OBJ = TTPCAStar.o TTPCHitPad.o TTPCLayout.o TTPCOrderedVolGroup.o TTPCPathVolume.o TTPCTRExPatAlgorithm.o TTPCTRExPatSubAlgorithm.o TTPCUnitVolume.o TTPCVolGroup.o TTPCVolGroupMan.o GasTPCDataLib.o TSimLoader.o GasTPCDataLibDict.o TTRExPattern.o TEventDisplay.o TTRExPatternDict.o TTRExHVCluster.o TTPCHelixPropagator.o TTPCLikelihoodMatch.o TTPCLikelihoodMerge.o TTPCLikFitPath.o TTPCQLikelihood.o TTPCSeeding.o TTPCTracking.o TTPCUtils.o TTrueHit.o TTRExPath.o TTRExJunction.o

ROOT_LIBS   = $(shell root-config --libs --glibs)
ROOT_CFLAGS = $(shell root-config --cflags)

#=========================================================
# Compile
#=========================================================

RunTRExRecon.exe: libpatrec.so RunTRExRecon.cxx
	$(CC) $(COPTS) $(ROOT_CFLAGS) $(ROOT_LIBS) $(INCLUDE) -I. -L. -lpatrec -lMinuit -o RunTRExRecon.exe RunTRExRecon.cxx

libpatrec.so: $(OBJ)
	$(CC) $(COPTS) -shared -o libpatrec.so $(OBJ)

#obj: $(OBJ)

%.o :	%.cxx
	$(CC) -c $(COPTS) $(ROOT_CFLAGS) $(INCLUDE) -I. -o $@ $<

GasTPCDataLibDict.cxx:
	rootcint GasTPCDataLibDict.cxx GasTPCDataLib.hxx LinkDef.hh

TTRExPatternDict.cxx: 
	rootcint TTRExPatternDict.cxx  TTRExPattern.hxx OutLink.hh

TTRExHVClusterDict.cxx:
	rootcint TTRExHVClusterDict.cxx TTRExHVCluster.hxx OutLink.hh

TTRExPathDict.cxx:
	rootcint TTRExPathDict.cxx TTRExPath.hxx OutLink.hh

clean:
	rm *.o *.so *.exe GasTPCDataLibDict.* TTRExPatternDict.* GasTPCDataLibDict_rdict.pcm TTRExPatternDict_rdict.pcm
