#include "TFile.h"
#include "TTree.h"
#include "TTRExPattern.hxx"

void DumpEvents(const char* fileName){
  TFile f(fileName);
  TTree* tr=(TTree*)f.Get("TPCRecon");
  trex::WritableEvent* ev=0;
  tr->SetBranchAddress("event",&ev);
  
  for(int i=0;i<tr->GetEntries();++i){
    tr->GetEntry(i);
    
    std::cout<<"Event "<<i<<":"<<std::endl;
    std::cout<<"contains "<<ev->patterns.size()<<" patterns"<<std::endl;

    for(std::vector<trex::WritablePattern>::iterator iPat=ev->patterns.begin();iPat!=ev->patterns.end();++iPat){
      std::cout<<"  Pattern containing "<<iPat->Paths.size()<<" paths and "<<iPat->Junctions.size()<<" junctions"<<std::endl;
      for(std::vector<std::vector<trex::TTPCHitPad> >::iterator iPath=iPat->Paths.begin();iPath!=iPat->Paths.end();++iPath){
	std::cout<<"    Path containing "<<iPath->size()<<" hits"<<std::endl;
      }
      for(std::vector<std::vector<trex::TTPCHitPad> >::iterator iJunct=iPat->Junctions.begin();iJunct!=iPat->Junctions.end();++iJunct){
	std::cout<<"    Junction containing "<<iJunct->size()<<" hits"<<std::endl;
      }
    }
  }
}
