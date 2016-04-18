#include "TSimLoader.hxx"

trex::TSimLoader::TSimLoader(std::string inputFile){

  fFile=new TFile(inputFile.c_str());

  fTree=(TTree*)fFile->Get("Tracking");

  fSimulDataBranch=0;
  fTree->SetBranchAddress("SimulData",&fSimulDataBranch);
}

void trex::TSimLoader::LoadEvent(unsigned int i){

  //Delete hits from last event
  for(std::vector<trex::TTPCHitPad*>::iterator hitPadIter=fHits.begin();hitPadIter!=fHits.end();++hitPadIter){
    delete *hitPadIter;
  }
  fHits.clear();

  fTree->GetEntry(i);

  HitCollection simHits = fSimulDataBranch->getTpcFidHits();

  for(HitCollection::iterator hitIter=simHits.begin();hitIter!=simHits.end();++hitIter){
    SDHit& hit=*hitIter;
    double eDep=hit.getEdep();
    TLorentzVector pos4=hit.getPosition();
    trex::TTPCHitPad* hitPadPtr=new TTPCHitPad(eDep,pos4);
    fHits.push_back(hitPadPtr);
  }
}

unsigned int trex::TSimLoader::GetNEvents(){
  return fTree->GetEntries();
}
