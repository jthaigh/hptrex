#include "TFile.h" 
#include "TTree.h"
#include "TTRExPattern.hxx"
#include "TTPCHitPad.hxx" 
#include "TH1.h"  


void PIDmacro(const char * pion, const char * proton){
  
  TFile * PionFile = new TFile(pion);
  TFile * ProtonFile = new TFile(proton);
  
  TTree * Pion_ReconTree = (TTree*)PionFile->Get("TPCRecon");
  TTree * Proton_ReconTree = (TTree*)ProtonFile->Get("TPCRecon");
  
  TH1D Pion_dEdx("Pi_dEdx", "Pi_dEdx", 100, 0, 0.00001);
  TH1D Pion_TrackLength("Pi_TrackLength", "Pi_TrackLength", 600, 0, 1200);
  TH1D Pion_ChargeSum("Pi_ChargeSum", "Pi_ChargeSum", 800, 0, 0.024);
  
  TH1D Proton_dEdx("P_dEdx", "P_dEdx", 100, 0, 0.00001); 
  TH1D Proton_TrackLength("P_TrackLength", "P_TrackLength", 600, 0, 1200);
  TH1D Proton_ChargeSum("P_ChargeSum", "P_ChargeSum", 800, 0, 0.024);
  
  int Pion_entries = Pion_ReconTree->GetEntries();
  int Proton_entries = Proton_ReconTree->GetEntries();
  
  trex::WritableEvent * Pion_event=0;
  trex::WritableEvent * Proton_event=0;
  
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

  TFile outf("SingeTrack_pi_p_PID.root", "RECREATE");

  Pion_dEdx.Write();
  Pion_TrackLength.Write();
  Pion_ChargeSum.Write();
  
  Proton_dEdx.Write();
  Proton_TrackLength.Write();
  Proton_ChargeSum.Write();

  Pion_dEdx.Draw();
  Proton_dEdx.Draw("same");
  
}

