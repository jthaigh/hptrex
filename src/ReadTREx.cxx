#include "TFile.h"
#include "TTree.h"
#include "TTRExPattern.hxx"
#include "TTPCHitPad.hxx"


void ReadTREx(const char * argv) {

  TFile * inFile = new TFile(argv);

  TTree * ReconTree = (TTree*)inFile->Get("TPCRecon");

  trex::TTRExEvent * event=0;

  ReconTree->SetBranchAddress("event", &event);

  for(int i=0;i<ReconTree->GetEntries();++i){
    ReconTree->GetEntry(i);
    std::vector<trex::TTRExPattern> pats = event->GetPatterns();
    int patSize = pats.size(); 
    std::cout << "Size of pattern vector: " << patSize << std::endl;

    for(int i=0; i<patSize; ++i){

    pats.at(i).Print();

    }

  }
}
