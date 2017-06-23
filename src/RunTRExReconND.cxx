#include "TTPCTRExPatAlgorithm.hxx"
#include "TTPCSeeding.hxx"
#include "TTPCTracking.hxx"
#include "TTPCLikelihoodMatch.hxx"
#include "TTPCLikelihoodMerge.hxx"
#include "TEventDisplay.hxx"
#include "TTPCHitPad.hxx"
#include "TSimLoaderND.hxx"
#include "TTrueHit.hxx"
#include "TTRExPattern.hxx"
#include "TTRExPIDAlgorithm.hxx"

#include <iostream>
#include <vector>
#include "TVector3.h"
#include "TFile.h"
#include "TClonesArray.h"


#include "TROOT.h"
#include "TRint.h"


int main(int argc, char** argv){
  
  gROOT->ProcessLine(".class trex::TTPCHitPad");

  trex::TSimLoaderND loader(argv[1]);

  const char * originalName = loader.GetFile()->GetName();
  std::cout << "You are processing File: " << originalName << std::endl; 
  string str(originalName);
  std::size_t found = str.rfind("/");
  str.erase(0,found+1);  
  const char * newName = str.c_str();
  char name[100];
  sprintf(name, "TRExRecon_%s", newName);  

  TFile * fOut = new TFile(name, "RECREATE"); 
  TTree * fReconTree = new TTree("TPCRecon", "TPCRecon");
  TTree * fTruthTree = new TTree("Truth", "Truth");
  fReconTree->SetDirectory(fOut);
  fTruthTree->SetDirectory(fOut);
  
  std::vector<trex::TTPCHitPad> * unused=0;
  trex::WritableEvent * outEvent=0;
  std::vector<trex::TTrueTrack> * outMCTracks=new std::vector<trex::TTrueTrack>;

  fReconTree->Branch("unusedHits", &unused);//, 64000, 1);
  fReconTree->Branch("event", &outEvent);//, 64000, 1);
  fTruthTree->Branch("tracks", &outMCTracks);//, 64000, 1);

  char merged_plots[100]; 
  char unmerged_plots[100];

  sprintf(merged_plots, "plots_merged_%s", newName);
  sprintf(unmerged_plots, "plots_unmerged_%s", newName);

  TFile fPlot(unmerged_plots,"RECREATE");
  TFile fPlotM(merged_plots, "RECREATE");

  //TFile fPlot("plots_unmerged.root","RECREATE");
  //TFile fPlotM("plots_merged.root","RECREATE");
  

  for(int i=0;i!= 100; ++i){ //loader.GetNEvents();++i){


    trex::TTRExEvent * event;
    
    loader.LoadEvent(i);
    
    std::vector<trex::TTPCHitPad*>& hitPads=loader.GetHits();
    std::vector<trex::TTPCHitPad*> usedHits;
    std::vector<TTrueHit*>& trueHits = loader.GetTrueHits();
    std::vector<trex::TTPCHitPad*> unusedTemp;
    
    outMCTracks->clear();
    std::vector<trex::TTrueTrack*>& evtTracks= loader.GetTrueTracks();

    for(auto iTrk=evtTracks.begin();iTrk!=evtTracks.end();++iTrk){
      outMCTracks->emplace_back(**iTrk);
      std::cout<<"Writing out an MC track"<<std::endl;
    }
    

    event = new trex::TTRExEvent();
        
    trex::TTPCTRExPatAlgorithm trexAlg;
    trex::TTPCSeeding seedingAlgo;
    trex::TTPCTracking trackingAlgo;
    trex::TTPCLikelihoodMatch matchAlgo;
    trex::TTPCLikelihoodMerge mergeAlgo;
    //Initialise PID Algo here
    trex::TTRExPIDAlgorithm pidAlgo;
    trex::TEventDisplay evDisp(&fPlot,i);
    trex::TEventDisplay evDispM(&fPlotM,i);

    trex::TTPCLayout layout;
    trexAlg.Process(hitPads,usedHits,unusedTemp,trueHits,event,layout); 
    
    for(auto iHit=unusedTemp.begin();iHit!=unusedTemp.end();++iHit){
      unused->emplace_back(**iHit);
    }

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
    std::cout <<"Running PID calculator..." <<std::endl;
    //Run PID Algo here
    pidAlgo.Process(mergedEvt->GetPatterns());

    evDisp.Process(hitPads,trueHits,event,layout);    
    evDispM.Process(hitPads,trueHits,mergedEvt,layout);    
    
    outEvent->FillFromEvent(*mergedEvt);
    fReconTree->Fill();
    fTruthTree->Fill();
    
    unused->clear();            
    delete event;
    delete mergedEvt;
    
  }
  
  fOut->Write();
  fOut->Close();
  fPlot.Write();
  fPlot.Close();
  fPlotM.Write();
  fPlotM.Close();
}

