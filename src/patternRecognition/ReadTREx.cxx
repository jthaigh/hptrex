#include "TFile.h"
#include "TTree.h"
#include "TTRExPattern.hxx"
#include "TTPCHitPad.hxx"


void ReadTREx(const char * argv) {

  TFile * inFile = new TFile(argv);

  TTree * ReconTree = (TTree*)inFile->Get("TPCRecon");

  trex::TTRExEvent * event;

  ReconTree->SetBranchAddress("event", &event);

  std::vector<trex::TTRExPattern> pats = event->GetPatterns();


  int patSize = pats.size(); 
  std::cout << "Size of pattern vector: " << patSize << std::endl;

  //for(int i=0; i<patSize; ++i){

  //pats.at(i).Print();

  //}

  //for(auto ipat=pats.begin(); ipat!=pats.end(); ++ipat){
 
  //std::vector<std::vector<trex::TTPCHitPad> > paths = ipat->GetPaths();
  //std::vector<std::vector<trex::TTPCHitPad> > juncts = ipat->GetJunctions();

  //}
}

