#include "TFile.h" 
#include "TTree.h"
#include "../TTRExPattern.hxx"
#include "../TTPCHitPad.hxx" 
#include "TH1.h"  
#include "TH2.h"
#include "TChain.h"
#include "TFileInfo.h"

void PIDmacro(){

  TChain * ReconTree = new TChain("TPCRecon");
  //ReconTree->Add("TRExRecon_plot_p5_vs_pi10_1_voxelsIdeal_234e-5m_234e-5m.root");
  //ReconTree->Add("");
  //ReconTree->Add("");
  //ReconTree->Add("");
  //ReconTree->Add("");

  std::cout << "We are reaching this 1" << std::endl;

  //declare histograms
  TH1D dEdx("dEdx", "dEdx", 50, 0, 0.00002);
  TH1D TrackLength("TrackLength", "TrackLength", 300, 0, 1200);
  TH1D ChargeSum("ChargeSum", "ChargeSum", 200, 0, 0.024);

  TH1D Cleanliness("Cleanliness", "Cleanliness", 100, 0, 1);
  
  TH1D* CleanlinessVsM[21];
  TH1D* CompletenessVsM[21];
  for(int i=0; i<21; ++i){
    CleanlinessVsM[i] = new TH1D("Cleanliness" , "Cleanliness", 100, 0, 1);   
    CompletenessVsM[i] = new TH1D("Completeness", "Completeness", 100, 0, 1);
  }
  
  TH1D Completeness("Completeness", "Completeness", 100, 0, 1);

  TH1D MeanCleanliness("MeanCleanliness", "MeanCleanliness", 20, 0, 20);
  TH1D MeanCompleteness("MeanCompleteness", "MeanCompleteness", 20, 0, 20);
  
  TH1D TrackEfficiency("TrackEfficiency", "TrackEfficiency", 20, 0, 20);
  TH1D ProtonEfficiency("ProtonEfficiency", "ProtonEfficicency",20, 0, 20); 
  TH1D ProtonPurity("ProtonPurity", "ProtonPurity", 20,0,20);
  
  TH2D Multiplicities("Multiplicities", "Multiplicities", 20, 0, 20, 20, 0, 20); 

  int entries = ReconTree->GetEntries();
  
  trex::WritableEvent * event=0;

  double pion_veto = 3.5e-6;

  ReconTree->SetBranchAddress("event", &event);
   

  std::cout << "check1" << std::endl;

  //loop over events
  
  int goodProtons[21];
  int badProtons[21];
  int missedProtons[21];

  int totalNumberOfPaths[21];
  int totalNumberOfGoodTracks[21];

  for(int i=0; i<entries; ++i) {
    
    ReconTree->GetEntry(i);

    std::cout << "check2" << std::endl;

    int recoMulti = event->RecoMultiplicity;
    int trueMulti = event->TrueMultiplicity;
    
    if (trueMulti >20){
      continue;
    }

    Multiplicities.Fill(trueMulti, recoMulti);
       
    std::vector<trex::WritablePattern> pats = event->patterns;
    for(int j=0; j<pats.size(); ++j){
      

      std::cout << "check3" << std::endl;

      std::vector<double> dEdx_vec = pats[j].dEdx;
      std::vector<double> TrackLength_vec = pats[j].TrackLength;
      std::vector<double> ChargeSum_vec = pats[j].ChargeSum;
      //declare more vectors
      std::vector<double> Cleanliness_vec = pats[j].TrackCleanliness;
      std::vector<double> Completeness_vec = pats[j].TrackCompleteness;
      std::vector<int> ProOrPi_vec = pats[j].ProOrPi;
      std::vector<int> PID_vec = pats[j].PID;

      for(int k=0; k<dEdx_vec.size(); ++k){
	
	std::cout << "check4" << std::endl;

	totalNumberOfPaths[trueMulti]+=1;
	
	dEdx.Fill(dEdx_vec[k]);
	TrackLength.Fill(TrackLength_vec[k]);
	ChargeSum.Fill(ChargeSum_vec[k]);
	Cleanliness.Fill(Cleanliness_vec[k]);
	Completeness.Fill(Completeness_vec[k]);
	
	if (Cleanliness_vec[k]>0.9 && Completeness_vec[k] > 0.8){
	  totalNumberOfGoodTracks[trueMulti]+=1;
	}
	
	CleanlinessVsM[trueMulti]->Fill(Cleanliness_vec[k]);
	CompletenessVsM[trueMulti]->Fill(Completeness_vec[k]); 
	
	if(ProOrPi_vec[k]==0 && PID_vec[k]==-1){
	  goodProtons[trueMulti]+=1;
	}else if(ProOrPi_vec[k]==1 && PID_vec[k]==-1){
	  badProtons[trueMulti]+=1;
	}else if(ProOrPi_vec[k]==0 && PID_vec[k]==1){
	  missedProtons[trueMulti]+=1;
	}
	
	std::cout << "check4.1" << std::endl;
	//fill more histograms
	
      }      

      std::cout << "check4.2" << std::endl;
    }

    std::cout << "check4.3" << std::endl;
  } 
  

  std::cout << "check5" << std::endl;

  for(int i=1; i<21; ++i) {
    
    std::cout << "check6" << std::endl;
    
    MeanCleanliness.SetBinContent(i,CleanlinessVsM[i]->GetMean());
    MeanCleanliness.SetBinError(i,CleanlinessVsM[i]->GetRMS());
    MeanCompleteness.SetBinContent(i, CompletenessVsM[i]->GetMean());
    MeanCompleteness.SetBinError(i, CompletenessVsM[i]->GetRMS());
    TrackEfficiency.SetBinContent(i, (double)totalNumberOfGoodTracks[i]/(double)totalNumberOfPaths[i]);
    
    int numberOfTrueProtons = goodProtons[i]+missedProtons[i];
    int numberOfRecoProtons = goodProtons[i]+badProtons[i];
    
    ProtonEfficiency.SetBinContent(i,(double)goodProtons[i]/(double)numberOfTrueProtons);
    ProtonPurity.SetBinContent(i, (double)goodProtons[i]/(double)numberOfRecoProtons);

    std::cout << "check7" << std::endl;

  }
  
  
  dEdx.SetLineColor(kBlue);
  TrackLength.SetLineColor(kBlue);
  ChargeSum.SetLineColor(kRed);
  
  std::cout << "check8" << std::endl;

  TFile outf("PID.root", "RECREATE");
  
  dEdx.Write();
  TrackLength.Write();
  ChargeSum.Write();
  Cleanliness.Write();
  Completeness.Write();
  MeanCleanliness.Write();
  MeanCompleteness.Write();
  TrackEfficiency.Write();
  ProtonEfficiency.Write();
  ProtonPurity.Write();
  Multiplicities.Write();
  //write more histograms


}

