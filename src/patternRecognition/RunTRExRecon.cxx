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

//#ifndef BOOST_ALL_DYN_LINK
//#   define BOOST_ALL_DYN_LINK
//#endif 

//#include <boost/filesystem.hpp>

int main(int argc, char** argv){
  
  gROOT->ProcessLine(".class trex::TTPCHitPad");

  trex::TSimLoader loader(argv[1]);

  TFile fOut("plots.root","RECREATE");
  //TFile fFile("recon.root", "RECREATE");
  //TFile fFile((TFile*)loader.GetFile()->Clone());
  //loader.GetFile()->Copy(fFile);
  const char * originalName = loader.GetFile()->GetName();
  std::cout << "You are processing File: " << originalName << std::endl; 
  TTree * fReconTree=(TTree*)loader.GetReconTree();
  
  std::vector<trex::TTPCHitPad> * unused;
  trex::TTRExEvent * event;

  fReconTree->Branch("unusedHits", &unused, 64000, 1);
  fReconTree->Branch("event", &event, 64000, 1);

  for(int i=0;i<500;++i){//loader.GetNEvents();++i){
    
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


    for(int i=0; i<unused->size();++i){

      std::cout << "Unused hits have content: " << std::endl;
      std::cout << "____________________________" << std::endl;
      unused->at(i).Print();
      std::cout << "____________________________" << std::endl;

    }
    
    std::cout << "Event " << i << " actually contains something: " << event->GetPatterns().size() << std::endl;
    
    
    fReconTree->Fill();
    delete event;
    delete unused;    
  }
  
  fReconTree->Print();
  

  //fOut.Write();
  //fOut.Close();

  //fFile.cd();


  //TFile treeFile("TreeFile.root", "RECREATE");
  //treeFile.cd();

  //fReconTree->Write();

  //treeFile.Write();
  //treeFile.Close();

  
  //std::set<char> delims{'/'};  
  //std::vector<std::string> path = splitpath(originalName, delims);
  //std::cout << path.back() << std::endl;
  
  
  //boost::filesystem::path p(originalName);
  //std::cout << "filename and extension : " << p.filename() << std::endl; // file.ext
  //std::cout << "filename only          : " << p.stem() << std::endl;     // file


  
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
  loader.GetFile()->Copy(*fFile);

  //fFile.SetName(name);
  std::cout << "This File has a Name: " << fFile->GetName() << std::endl;
 
  fReconTree->Write();
  fFile->Write(name);
  fFile->Close();
}
    
