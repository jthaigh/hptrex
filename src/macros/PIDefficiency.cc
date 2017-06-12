#include "TFile.h" 
#include "TTree.h"
#include "TTRExPattern.hxx"
#include "TTPCHitPad.hxx" 
#include "TH1.h"  
#include "TChain.h"
#include "TFileInfo.h"

void PIDmacro(const char * pion, const char * proton){
  
  //TChain * PionFile = new TChain("TPCRecon");
  //PionFile->Add(pion);

  //TChain * ProtonFile = new TChain("TPCRecon");
  //ProtonFile->Add(proton); 
  
  TFile * PionFile = new TFile(pion);
  TFile * ProtonFile = new TFile(proton);
    
  TTree * Pion_ReconTree = (TTree*)PionFile->Get("TPCRecon");
  TTree * Proton_ReconTree = (TTree*)ProtonFile->Get("TPCRecon");

  std::cout << "We are reaching this 1" << std::endl;

  //TTree * Pion_ReconTree = (TTree*)PionFile->GetTree();
  //TTree * Proton_ReconTree = (TTree*)ProtonFile->GetTree(); 

  std::cout << "We are able to get the Tree" << std::endl;

  TH1D Pion_dEdx("Pi_dEdx", "Pi_dEdx", 50, 0, 0.00002);
  TH1D Pion_TrackLength("Pi_TrackLength", "Pi_TrackLength", 300, 0, 1200);
  TH1D Pion_ChargeSum("Pi_ChargeSum", "Pi_ChargeSum", 200, 0, 0.024);
  
  TH1D Proton_dEdx("P_dEdx", "P_dEdx", 50, 0, 0.00002); 
  TH1D Proton_TrackLength("P_TrackLength", "P_TrackLength", 300, 0, 1200);
  TH1D Proton_ChargeSum("P_ChargeSum", "P_ChargeSum", 200, 0, 0.024);
  
  int Pion_entries = Pion_ReconTree->GetEntries();
  int Proton_entries = Proton_ReconTree->GetEntries();
  
  trex::WritableEvent * Pion_event=0;
  trex::WritableEvent * Proton_event=0;

  double pion_veto = 3.5e-6;

  int protons_passed = 0;
  int pions_passed = 0; 

  int protons_total = 0;
  int pions_total = 0;

  
  Pion_ReconTree->SetBranchAddress("event", &Pion_event);
  Proton_ReconTree->SetBranchAddress("event", &Proton_event);
  
  //loop over pion events
  
  for(int i=0; i<Pion_entries; ++i) {
    
    Pion_ReconTree->GetEntry(i);
    std::vector<trex::WritablePattern> pats = Pion_event->patterns;
    for(int j=0; j<pats.size(); ++j){
      
      std::vector<double> dEdx_vec = pats[j].dEdx;
      std::vector<double> TrackLength_vec = pats[j].TrackLength;
      std::vector<double> ChargeSum_vec = pats[j].ChargeSum;
      
      for(int k=0; k<dEdx_vec.size(); ++k){

	if(dEdx_vec[k] > pion_veto){pions_passed += 1;}
	pions_total += 1;
	Pion_dEdx.Fill(dEdx_vec[k]);
	Pion_TrackLength.Fill(TrackLength_vec[k]);
	Pion_ChargeSum.Fill(ChargeSum_vec[k]);
      }      
    } 
  }
  
  for(int i=0; i<Proton_entries; ++i) {
    Proton_ReconTree->GetEntry(i);
    std::vector<trex::WritablePattern> pats = Proton_event->patterns;
    for(int j=0; j<pats.size(); ++j){
    
      std::vector<double> dEdx_vec = pats[j].dEdx;
      std::vector<double> TrackLength_vec = pats[j].TrackLength;
      std::vector<double> ChargeSum_vec = pats[j].ChargeSum;

      for(int k=0; k<dEdx_vec.size(); ++k){
	
	if(dEdx_vec[k] > pion_veto){protons_passed += 1;}
	protons_total += 1;
	Proton_dEdx.Fill(dEdx_vec[k]);
	Proton_TrackLength.Fill(TrackLength_vec[k]);
	Proton_ChargeSum.Fill(ChargeSum_vec[k]);
      }
    }
  }
  
  
  Pion_dEdx.SetLineColor(kBlue);
  Pion_TrackLength.SetLineColor(kBlue);
  Pion_ChargeSum.SetLineColor(kBlue);
  
  Proton_dEdx.SetLineColor(kRed);
  Proton_TrackLength.SetLineColor(kRed);
  Proton_ChargeSum.SetLineColor(kRed);

  TFile outf("pi_p_PID.root", "RECREATE");

  Pion_dEdx.Write();
  Pion_TrackLength.Write();
  Pion_ChargeSum.Write();
  
  Proton_dEdx.Write();
  Proton_TrackLength.Write();
  Proton_ChargeSum.Write();

  Pion_dEdx.Draw();
  Proton_dEdx.Draw("same");
  
  std::cout << "TOTAL PIONS: " << pions_total << " OF WHICH PASSED " << pions_passed << std::endl;

  std::cout << "TOTAL PROTONS: " << protons_total << " OF WHICH PASSED " << protons_passed << std::endl;

  std::cout << "PROTON PURITY: " << pions_passed/protons_passed << std::endl; 
  std::cout << "PROTON EFFICIENCY: " << protons_passed/protons_total << std::endl;
}

