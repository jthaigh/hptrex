
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


  //make pointers for output variables to be filled by Process()
  std::vector<trex::TTPCHitPad> * unused;
  trex::TTRExEvent * event;

  fReconTree->Branch("unusedHits", &unused, 64000, 1);
  fReconTree->Branch("event", &event, 64000, 1);

  for(int i=0;i<2000;++i){//loader.GetNEvents();++i){
    
    loader.LoadEvent(i);
    
    std::vector<trex::TTPCHitPad*>& hitPads=loader.GetHits();
    std::vector<trex::TTPCHitPad*> usedHits;
    std::vector<trex::TTPCHitPad*> unusedHits;
    std::vector<TTrueHit*>& trueHits = loader.GetTrueHits();

    event = new trex::TTRExEvent();
    unused = new std::vector<trex::TTPCHitPad>();
    
    std::cout << "True hits contains: " << trueHits.size() << " entries. "<< std::endl;

    //loader.DrawDetector();
    
    trex::TTPCTRExPatAlgorithm trexAlg(&fOut);
    std::cout<<"EVERYTHING LOADED! - NOW ATTEMPTING TO PROCESS"<<std::endl;
    trexAlg.Process(hitPads,usedHits,unused,trueHits,event); 
          
    //}
    
    std::cout << "SIZE OF UNUSED HITS VECTOR " << unused->size() << std::endl;

    
    for(int j=0; j<unused->size();++j){
      
      std::cout << "Unused hits have content: " << std::endl;
      std::cout << "____________________________" << std::endl;
      unused->at(j).Print();
      std::cout << "____________________________" << std::endl;

    }
    
    std::cout << "Event " << i << " contains  "<< event->GetPatterns().size() << " Patterns" <<std::endl;
    
   
    for(int k=0; k<event->GetPatterns().size(); ++k){
      
      std::cout << "PATTERN NUMBER " << k << std::endl;

      //Print out Path information for debugging
      for(int l=0; l<event->GetPatterns().at(k).GetPaths().size(); ++l){	
	std::cout << "PATH " << l << " contains the following: " << std::endl;
	std::cout << "____________________________" << std::endl;
	event->GetPatterns().at(k).GetPaths().at(l).Print();
	std::cout << "____________________________" << std::endl;
      } 
      
      //Print out Junction Information for debugging
      for(int l=0; l<event->GetPatterns().at(k).GetJunctions().size(); ++l){
	
	std::cout << "JUNCTION " << l << "contains the following: " << std::endl;
	std::cout << "____________________________" << std::endl;
	event->GetPatterns().at(k).GetJunctions().at(l).Print();
	std::cout << "____________________________" << std::endl;
	
      }
      
      
      std::cout << "CONNECTED my Map of size " << event->GetPatterns().at(k).GetMap().size() << std::endl;
      std::cout << "Map has following entries: " << std::endl;
      std::cout << "____________________________" << std::endl;

      //Print out Map information for debugging
      for (int l=0; l<event->GetPatterns().at(k).GetMap().size(); ++l){
	
	for (int m=0; m<event->GetPatterns().at(k).GetMap().at(l).size(); ++m){
	  std::cout << event->GetPatterns().at(k).GetMap().at(l).at(m)<< std::endl;	
	}
      }      
      std::cout << "____________________________" << std::endl;
    }

    


    //DO TRACKING HERE using event object filled above
    





    //fReconTree->Fill();
    
    delete event;
    delete unused;        
    
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

