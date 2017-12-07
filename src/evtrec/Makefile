#=========================================================
# Makefile
#=========================================================

CC	= g++
COPTS	= -Wall -g -fPIC -DLINUX -O2 -std=c++11
# -Wl,-rpath,$(shell brew --prefix)/lib

OBJ = EventRecord.o MCGenInfo.o MCHit.o MCParticle.o MCTrack.o RecoParticle.o RootFileReader.o RootFileWriter.o Vector3D.o Vector4D.o EventRecordDict.o

ROOT_LIBS   = $(shell root-config --libs --glibs)
ROOT_CFLAGS = $(shell root-config --cflags)

#=========================================================
# Compile
#=========================================================

libevtrecord.so: $(OBJ)
	$(CC) $(COPTS) -shared -o libevtrecord.so $(OBJ)

%.o :	%.cxx
	$(CC) -c $(COPTS) $(ROOT_CFLAGS) $(INCLUDE) -I. -o $@ $<

EventRecordDict.cxx:
	rootcint EventRecordDict.cxx EventRecord.h MCGenInfo.h MCParticle.h MCTrack.h MCHit.h RecoParticle.h Vector3D.h Vector4D.h LinkDef.h

clean:
	rm -f *.o *.so EventRecordDict.*
