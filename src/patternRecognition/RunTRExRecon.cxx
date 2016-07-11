#include "TTPCTRExPatAlgorithm.hxx"
#include "TTPCHitPad.hxx"
#include "TSimLoader.hxx"

#include <iostream>
#include <vector>
#include "TVector3.h"
#include "TFile.h"

int main(int argc,const char** argv){
  
  trex::TSimLoader loader(argv[1]);

  TFile fOut("plots.root","RECREATE");

  for(int i=0;i<1000;++i){//loader.GetNEvents();++i){
    
    loader.LoadEvent(i);
    std::vector<trex::TTPCHitPad*>& hitPads=loader.GetHits();
    std::vector<trex::TTPCHitPad*> usedHits;
    std::vector<trex::TTPCHitPad*> unusedHits;

    trex::TTPCTRExPatAlgorithm trexAlg(&fOut);
    std::cout<<"EVERYTHING LOADED! - NOW ATTEMPTING TO PROCESS"<<std::endl;
    trexAlg.Process(hitPads,usedHits,unusedHits);
  }

  fOut.Write();
}
    
