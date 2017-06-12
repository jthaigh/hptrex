#include "TFile.h" 
#include "TTree.h" 
#include "TH1.h" 

void AddPIDhists(const char * first, const char * second){

  TFile * firstFile = new TFile(first);
  TFile * secondFile = new TFile(second);

  TH1D * first_Pi_dEdx = (TH1D*)firstFile->Get("Pi_dEdx");
  TH1D * first_Pi_TrackLength = (TH1D*)firstFile->Get("Pi_TrackLength");
  TH1D * first_Pi_ChargeSum = (TH1D*)firstFile->Get("Pi_TrackLength");

  TH1D * second_Pi_dEdx = (TH1D*)secondFile->Get("Pi_dEdx");
  TH1D * second_Pi_TrackLength = (TH1D*)secondFile->Get("Pi_TrackLength");
  TH1D * second_Pi_ChargeSum = (TH1D*)secondFile->Get("Pi_TrackLength");

  TH1D * first_P_dEdx = (TH1D*)firstFile->Get("P_dEdx"); 
  TH1D * first_P_TrackLength = (TH1D*)firstFile->Get("P_TrackLength");
  TH1D * first_P_ChargeSum = (TH1D*)firstFile->Get("P_TrackLength");

  TH1D * second_P_dEdx = (TH1D*)secondFile->Get("P_dEdx");
  TH1D * second_P_TrackLength = (TH1D*)secondFile->Get("P_TrackLength");
  TH1D * second_P_ChargeSum = (TH1D*)secondFile->Get("P_TrackLength");



  TH1D Pion_dEdx("Pi_dEdx", "Pi_dEdx", 50, 0, 0.00002);
  TH1D Pion_TrackLength("Pi_TrackLength", "Pi_TrackLength", 300, 0, 1200);
  TH1D Pion_ChargeSum("Pi_ChargeSum", "Pi_ChargeSum", 200, 0, 0.024); 


  Pion_dEdx.Add(first_Pi_dEdx, second_Pi_dEdx);
  Pion_TrackLength.Add(first_Pi_TrackLength, second_Pi_TrackLength);
  Pion_ChargeSum.Add(first_Pi_ChargeSum, second_Pi_ChargeSum);

  TH1D Proton_dEdx("P_dEdx", "P_dEdx", 50, 0, 0.00002);
  TH1D Proton_TrackLength("P_TrackLength", "P_TrackLength", 300, 0, 1200);
  TH1D Proton_ChargeSum("P_ChargeSum", "P_ChargeSum", 200, 0, 0.024);

  Proton_dEdx.Add(first_P_dEdx, second_P_dEdx);
  Proton_TrackLength.Add(first_P_TrackLength, second_P_TrackLength);
  Proton_ChargeSum.Add(first_P_ChargeSum, second_P_ChargeSum);


  Pion_dEdx.SetLineColor(kBlue);
  Pion_TrackLength.SetLineColor(kBlue);
  Pion_ChargeSum.SetLineColor(kBlue);

  Proton_dEdx.SetLineColor(kRed);
  Proton_TrackLength.SetLineColor(kRed);
  Proton_ChargeSum.SetLineColor(kRed);

  TFile outf("ADDED_UP_pi_p_PID.root", "RECREATE");

  Pion_dEdx.Write();
  Pion_TrackLength.Write();
  Pion_ChargeSum.Write();

  Proton_dEdx.Write();
  Proton_TrackLength.Write();
  Proton_ChargeSum.Write();


}
