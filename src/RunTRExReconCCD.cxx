#include "TTPCTRExPatAlgorithm.hxx"
#include "TTPCLinearMerge.hxx"
#include "TEventDisplay.hxx"
#include "TTPCHitPad.hxx"
#include "TSimLoaderCCD.hxx"
#include "TTrueHit.hxx"
#include "TTRExPattern.hxx"
#include "TTRExPIDAlgorithm.hxx"
#include "TTrueTrack.hxx"

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
  fReconTree->SetDirectory(fOut);
  
  std::vector<trex::TTPCHitPad> * unused=0;
  trex::WritableEvent * outEvent=0;

  fReconTree->Branch("unusedHits", &unused);
  fReconTree->Branch("event", &outEvent);

  char merged_plots[100];
  char unmerged_plots[100];

  sprintf(merged_plots, "plots_merged_%s", newName);
  sprintf(unmerged_plots, "plots_unmerged_%s", newName);

  TFile fPlot(unmerged_plots,"RECREATE");
  TFile fPlotM(merged_plots, "RECREATE");
    
  //TFile fPlot("plots_unmerged.root","RECREATE");
  //TFile fPlotM("plots_merged.root","RECREATE");
  
  for(int i=0;i<loader.GetNEvents();++i){

    trex::TTRExEvent * event;
    
    loader.LoadEvent(i);
    
    std::vector<trex::TTPCHitPad*>& hitPads=loader.GetHits();
    std::vector<trex::TTPCHitPad*> usedHits;
    std::vector<TTrueHit*>& trueHits = loader.GetTrueHits();
    std::vector<trex::TTPCHitPad*> unusedTemp;

    int TrueMultiplicity=loader.GetTrueMultiplicity();

    event = new trex::TTRExEvent();
        
    trex::TTPCTRExPatAlgorithm trexAlg;
    trex::TTPCLinearMerge mergeAlgo;
    //Initialise PID Algo here
    trex::TTRExPIDAlgorithm pidAlgo;
    std::cout << "Now setting up event displays" << std::endl;
    trex::TEventDisplay evDisp(&fPlot,i);
    trex::TEventDisplay evDispM(&fPlotM,i);

    trex::TTPCLayout layout;

    trexAlg.Process(hitPads,usedHits,unusedTemp,trueHits,event,layout); 
    

    
    for(auto iHit=unusedTemp.begin();iHit!=unusedTemp.end();++iHit){
      unused->emplace_back(**iHit);
    }

    std::cout<<"Pattern recognition produced "<<event->GetPatterns().size()<<" patterns"<<std::endl;

    trex::TTRExEvent* mergedEvt=new trex::TTRExEvent;
    std::cout<<"Running merging..."<<std::endl;
    mergeAlgo.Process(event->GetPatterns(),mergedEvt->GetPatterns());
    
    std::cout <<"Running PID calculator..." <<std::endl;

    //Run PID Algo here
    pidAlgo.Process(mergedEvt->GetPatterns());
    
    std::cout << "Now filling Event Displays" << std::endl;
    evDisp.Process(hitPads,trueHits,event,layout);    
    evDispM.Process(hitPads,trueHits,mergedEvt,layout);    
    
    std::cout << "Now building output Event for i/o" << std::endl;
    outEvent->FillFromEvent(*mergedEvt);
    outEvent->SetTrueMultiplicity(TrueMultiplicity);

    std::cout << "Now filling Recon Tree for output" << std::endl;
    fReconTree->Fill();
    
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

