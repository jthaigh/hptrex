
#include "TTPCTRExPatAlgorithm.hxx"
#include "TTPCSeeding.hxx"
#include "TTPCTracking.hxx"
#include "TTPCLikelihoodMatch.hxx"
#include "TTPCLikelihoodMerge.hxx"
#include "TEventDisplay.hxx"
#include "TTPCHitPad.hxx"
#include "TSimLoader.hxx"
#include "TTrueHit.hxx"
#include "TTRExPattern.hxx"

#include <iostream>
#include <vector>
#include "TVector3.h"
#include "TFile.h"
#include "TClonesArray.h"


#include "TROOT.h"
#include "TRint.h"


int main(int argc, char** argv){
  
  gROOT->ProcessLine(".class trex::TTPCHitPad");

  trex::TSimLoader loader(argv[1]);

  TFile fOut("plots.root","RECREATE");
  
  const char * originalName = loader.GetFile()->GetName();
  std::cout << "You are processing File: " << originalName << std::endl; 
  
  string str(originalName);
  std::size_t found = str.rfind("/");

  str.erase(0,found+1);
  
  //  std::cout << "This is what found says: " << found << std::endl;
  //std::cout << "This ist what str says: " << str << std::endl;
  
  const char * newName = str.c_str();
  
  char name[100];
  sprintf(name, "TRExRecon_%s", newName);
  std:: cout << "New file will be called: " << name << std::endl;
  
  TFile * fFile = new TFile(name, "RECREATE"); 

  //string inputFile = argv[1];
  //TFile * fFile = new TFile(inputFile.c_str(), "UPDATE");
                                                                                                                        
  //std::cout << "This File has a Name: " << fFile->GetName() << std::endl;


  TTree * fReconTree = new TTree("TPCRecon", "TPCRecon");
  fReconTree->SetDirectory(fFile);
  //TTree * fReconTree=(TTree*)loader.GetReconTree();


  //make pointers for output variables to be filled by Process()
  std::vector<trex::TTPCHitPad*> unused;
  trex::TTRExEvent * event;

  fReconTree->Branch("unusedHits", &unused, 64000, 1);
  fReconTree->Branch("event", &event, 64000, 1);

  for(int i=0;i<loader.GetNEvents();++i){
    
    loader.LoadEvent(i);
    
    std::vector<trex::TTPCHitPad*>& hitPads=loader.GetHits();
    std::vector<trex::TTPCHitPad*> usedHits;
    std::vector<TTrueHit*>& trueHits = loader.GetTrueHits();

    event = new trex::TTRExEvent();
    
    //std::cout << "True hits contains: " << trueHits.size() << " entries. "<< std::endl;

    //loader.DrawDetector();
    
    trex::TTPCTRExPatAlgorithm trexAlg;
    trex::TTPCSeeding seedingAlgo;
    trex::TTPCTracking trackingAlgo;
    trex::TTPCLikelihoodMatch matchAlgo;
    trex::TTPCLikelihoodMerge mergeAlgo;
    trex::TEventDisplay evDisp(&fOut);

    trex::TTPCLayout layout;
    //    std::cout<<"EVERYTHING LOADED! - NOW ATTEMPTING TO PROCESS"<<std::endl;
    trexAlg.Process(hitPads,usedHits,unused,trueHits,event,layout); 

    //    if ( !event->GetPatterns().size())
    //  {
    //	continue;
    //  }

    for (auto pattern = event->GetPatterns().begin(); pattern != event->GetPatterns().end(); pattern++) {
      seedingAlgo.Process(*pattern);
      trackingAlgo.Process(*pattern);
    }

    trex::TTRExEvent* mergedEvt=new trex::TTRExEvent;
    
    matchAlgo.Process(event->GetPatterns());
    mergeAlgo.Process(event->GetPatterns(),mergedEvt->GetPatterns());

    evDisp.Process(hitPads,trueHits,event,layout);    

    //for(auto i=event->GetPatterns().begin();i!=event->GetPatterns().end();++i){
    //  i->Print();
    //}





    //fReconTree->Fill();
    
    delete event;
    //delete unused;        
    delete mergedEvt;
    
  }
  //fReconTree->Print();
  
  
  //fOut.Write();
  //fOut.Close();
  
  //fFile->cd();
  
  //std::cout << "Writing Tree to File" << std::endl;
  //fReconTree->Write();
  
  //fFile->Write();
  //fFile->Close();
}

