#include "TTPCTRExPatAlgorithm.hxx"
#include "TTPCHitPad.hxx"
#include "TSimLoader.hxx"

#include <iostream>
#include <vector>
#include "TVector3.h"

int main(int argc,const char** argv){
  
  trex::TSimLoader loader(argv[1]);

  for(int i=0;i<loader.GetNEvents();++i){
    
    loader.LoadEvent(i);
    std::vector<trex::TTPCHitPad*>& hitPads=loader.GetHits();
    std::vector<trex::TTPCHitPad*> usedHits;
    std::vector<trex::TTPCHitPad*> unusedHits;

    for(int i=0;i<hitPads.size();++i){
      TVector3 pos=hitPads[i]->GetPosition();
      std::cout<<"Hitpos: "<<pos.X()<<", "<<pos.Y()<<", "<<pos.Z()<<std::endl;
    }

    trex::TTPCTRExPatAlgorithm trexAlg;
    std::cout<<"1"<<std::endl;
    trexAlg.Process(hitPads,usedHits,unusedHits);
  }
}
    
