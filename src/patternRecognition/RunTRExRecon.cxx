#include "TTPCTRExPatAlgorithm.hxx"
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
  
  std::cout << "This is what found says: " << found << std::endl;
  std::cout << "This ist what str says: " << str << std::endl;
  
  const char * newName = str.c_str();
  
  char name[100];
  sprintf(name, "TRExRecon_%s", newName);
  std:: cout << "New file will be called: " << name << std::endl;
  
  TFile * fFile = new TFile(name, "RECREATE"); 

  //string inputFile = argv[1];
  //TFile * fFile = new TFile(inputFile.c_str(), "UPDATE");
                                                                                                                        
  std::cout << "This File has a Name: " << fFile->GetName() << std::endl;


  TTree * fReconTree = new TTree("TPCRecon", "TPCRecon");
  fReconTree->SetDirectory(fFile);
  //TTree * fReconTree=(TTree*)loader.GetReconTree();
  
  std::vector<trex::TTPCHitPad> * unused;
  trex::WritableEvent * outEvent=0;

  fReconTree->Branch("unusedHits", &unused);//, 64000, 1);
  fReconTree->Branch("event", &outEvent);//, 64000, 1);

  for(int i=0;i<2000;++i){//loader.GetNEvents();++i){
    
    loader.LoadEvent(i);
    
    std::vector<trex::TTPCHitPad*>& hitPads=loader.GetHits();
    std::vector<trex::TTPCHitPad*> usedHits;
    std::vector<trex::TTPCHitPad*> unusedHits;
    std::vector<TTrueHit*>& trueHits = loader.GetTrueHits();

    trex::TTRExEvent event;
    unused = new std::vector<trex::TTPCHitPad>();
    
    std::cout << "True hits contains: " << trueHits.size() << " entries. "<< std::endl;

    //loader.DrawDetector();
    
    trex::TTPCTRExPatAlgorithm trexAlg(&fOut);
    std::cout<<"EVERYTHING LOADED! - NOW ATTEMPTING TO PROCESS"<<std::endl;
    trexAlg.Process(hitPads,usedHits,unused,trueHits,&event); 
          
    //}
    
    std::cout << "SIZE OF UNUSED HITS VECTOR " << unused->size() << std::endl;

    
    for(int j=0; j<unused->size();++j){
      
      std::cout << "Unused hits have content: " << std::endl;
      std::cout << "____________________________" << std::endl;
      unused->at(j).Print();
      std::cout << "____________________________" << std::endl;

    }
    
    std::cout << "Event " << i << " actually contains something: " << event.GetPatterns().size() << std::endl;

    if(event.GetPatterns().size()!=0){std::cout << "EVENT DOES CONTAIN SOMETHING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << std::endl << std::endl;}
    
   
    for(int k=0; k<event.GetPatterns().size(); ++k){
      
      for(int l=0; l<event.GetPatterns().at(k).GetPaths().size(); ++l){
	
	for (int m=0; m<event.GetPatterns().at(k).GetPaths().at(l).size(); ++m){
	  
	  std::cout << "PATH hits have content: " << std::endl;
	  std::cout << "____________________________" << std::endl;
	  event.GetPatterns().at(k).GetPaths().at(l).at(m)->Print();
	  std::cout << "____________________________" << std::endl;
	  
	}
      } 
    }
    
    
    for(int k=0; k<event.GetPatterns().size(); ++k){
      
      for(int l=0; l<event.GetPatterns().at(k).GetJunctions().size(); ++l){
	
	for (int m=0; m<event.GetPatterns().at(k).GetJunctions().at(l).size(); ++m){
	  
	  std::cout << "JUNCTION hits have content: " << std::endl;
	  std::cout << "____________________________" << std::endl;
	  event.GetPatterns().at(k).GetJunctions().at(l).at(m)->Print();
	  std::cout << "____________________________" << std::endl;
	}
      }
    }
    
    

    /*    if(event->GetPatterns().size()==0){
    delete event;
    event=NULL;
    std::cout << "Event got deleted" << std::endl;
    }
    
    
    if(unused->size()==0){
    delete unused;
    unused=NULL;
    std::cout << "unused Hits got deleted" << std::endl;
    }*/
    
    outEvent->FillFromEvent(event);
    fReconTree->Fill();
    
    delete unused;        
    
  }
  fReconTree->Print();
  
  
  //fOut.Write();
  //fOut.Close();
  
  //fFile->cd();
  
  std::cout << "Writing Tree to File" << std::endl;
  fReconTree->Write();
  
  fFile->Write();
  fFile->Close();
}

