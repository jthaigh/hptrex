#ifndef TSIMLOADER_HXX
#define TSIMLOADER_HXX

//C++
#include <vector>

// ROOT
#include <TVector3.h>
#include "TFile.h"
#include "TTree.h"

//TREx
#include "TTPCHitPad.hxx"

//Sim data
#include "GasTPCDataLib.hxx"

namespace trex{

  class TSimLoader{

  public:

    TSimLoader(std::string inputFile);

    void LoadEvent(unsigned int i);

    unsigned int GetNEvents();

    inline std::vector<trex::TTPCHitPad*>& GetHits(){return fHits;}

  private:

    TFile* fFile;

    TTree* fTree;

    std::vector<trex::TTPCHitPad*> fHits;
    
    SimulData* fSimulDataBranch;

  };
}

#endif
