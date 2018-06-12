#include "TFile.h"
#include "TTree.h"
#include "../TTRExPattern.hxx"
#include "../TTPCHitPad.hxx"
#include "TH1.h"


void ReadTREx(const char * argv) {

  TFile * inFile = new TFile(argv);

  TTree * ReconTree = (TTree*)inFile->Get("TPCRecon");


  
  TH1D xPaths("xPosition", "xPosition", 1200, -600, 600);
  TH1D yPaths("yPosition", "yPosition", 1200, -600, 600);
  TH1D zPaths("zPosition", "zPosition", 1200, 0, 1200);

  
  TH1D xJuncts("xPosition", "xPosition", 1200, -600, 600);
  TH1D yJuncts("yPosition", "yPosition", 1200, -600, 600);
  TH1D zJuncts("zPosition", "zPosition", 1200, 0, 1200);
  
  //int entries = ReconTree->GetEntries();
  
  
  trex::TTRExEvent * event=0;
      
  ReconTree->SetBranchAddress("event", &event);
  
  for(int i=0;i<ReconTree->GetEntries();++i){
    ReconTree->GetEntry(i);
    std::vector<trex::TTRExPattern> pats = event->GetPatterns();
    int patSize = pats.size(); 
    std::cout << "Size of pattern vector: " << patSize << std::endl;
    
    for(int i=0; i<patSize; ++i){
      
      pats.at(i).Print();
      
      for (int j=0; j<pats.at(i).GetPaths().size(); ++j){
	
	for (int k=0; k<pats.at(i).GetPaths().at(j).size(); ++k){
	  
	  TVector3 pos = pats.at(i).GetPaths().at(j).at(k).GetPosition();
	  xPaths.Fill(pos.X());
	  yPaths.Fill(pos.Y());
	  zPaths.Fill(pos.Z());
	}
      }
      
      
      for (int l=0; l<pats.at(i).GetJunctions().size(); ++l){
	
	for ( int k=0; k<pats.at(i).GetJunctions().at(l).size(); ++k){
	  
	  TVector3 pos = pats.at(i).GetJunctions().at(l).at(k).GetPosition(); 
	  xJuncts.Fill(pos.X());
	  yJuncts.Fill(pos.Y());
	  zJuncts.Fill(pos.Z());
	}
	
      }
      
    }
    
  }

  TFile out("PositionDists.root", "RECREATE");
  
  xPaths.Write();
  yPaths.Write();
  zPaths.Write();
  
  xJuncts.Write();
  yJuncts.Write();
  zJuncts.Write();
  
}
