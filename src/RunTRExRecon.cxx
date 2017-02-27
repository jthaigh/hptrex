
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

  TFile fOut("plots_unmerged.root","RECREATE");
  TFile fOutM("plots_merged.root","RECREATE");
  
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
  
  //for now:
  //TFile * fFile = new TFile(name, "RECREATE"); 

  //string inputFile = argv[1];
  //TFile * fFile = new TFile(inputFile.c_str(), "UPDATE");
                                                                                                                        
  //std::cout << "This File has a Name: " << fFile->GetName() << std::endl;


  TTree * fReconTree = new TTree("TPCRecon", "TPCRecon");
  //fReconTree->SetDirectory(fFile);
  //TTree * fReconTree=(TTree*)loader.GetReconTree();


  //make pointers for output variables to be filled by Process()
  std::vector<trex::TTPCHitPad*> unused;
  trex::TTRExEvent * event;

  fReconTree->Branch("unusedHits", &unused, 64000, 1);
  fReconTree->Branch("event", &event, 64000, 1);

  for(int i=0;i!=loader.GetNEvents();++i){
    
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
    trex::TEventDisplay evDisp(&fOut,i);
    trex::TEventDisplay evDispM(&fOutM,i);

    trex::TTPCLayout layout;
    //    std::cout<<"EVERYTHING LOADED! - NOW ATTEMPTING TO PROCESS"<<std::endl;
    trexAlg.Process(hitPads,usedHits,unused,trueHits,event,layout); 

    //if ( !event->GetPatterns().size())
    //{
    //	continue;
    //}

    std::cout<<"Pattern recognition produced "<<event->GetPatterns().size()<<" patterns"<<std::endl;

    for (auto pattern = event->GetPatterns().begin(); pattern != event->GetPatterns().end(); pattern++) {
      std::cout<<"Running seeding on a pattern"<<std::endl;
      seedingAlgo.Process(*pattern);
      std::cout<<"Running tracking on a pattern"<<std::endl;
      trackingAlgo.Process(*pattern);
    }

    trex::TTRExEvent* mergedEvt=new trex::TTRExEvent;
    std::cout<<"Running matching..."<<std::endl;
    matchAlgo.Process(event->GetPatterns());
    std::cout<<"Running merging..."<<std::endl;
    mergeAlgo.Process(event->GetPatterns(),mergedEvt->GetPatterns());

    evDisp.Process(hitPads,trueHits,event,layout);    
    evDispM.Process(hitPads,trueHits,mergedEvt,layout);    

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

