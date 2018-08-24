#include "TFile.h" 
#include "TTree.h"
#include "../TTRExPattern.hxx"
#include "../TTPCHitPad.hxx" 
#include "TH1.h"  
#include "TH2.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TFileInfo.h"
#include "math.h"
#include "TStyle.h"


void PIDmacro(const char* inputFile){

  TChain * ReconTree = new TChain("TPCRecon");
 
  

/*  
  ReconTree->Add("~/TREx/hptrex/src/TRExRecon_plot_pBiasedElasticProton7.7_vs_piUnbiasedPion19.3_voxelsIdeal_234e-5m_234e-5m.root");
  ReconTree->Add("~/TREx/hptrex/src/TRExRecon_plot_pBiasedInelasticProton7.7_vs_piUnbiasedPion19.3_voxelsIdeal_234e-5m_234e-5m.root");
  ReconTree->Add("~/TREx/hptrex/src/TRExRecon_plot_pUnbiasedProton7.7_vs_piUnbiasedPion19.3_voxelsIdeal_234e-5m_234e-5m.root");
*/

  //ReconTree->Add("~/TREx/hptrex/src/Data/tobyData/MergedFiles280318/Processed/277_554/TRExRecon_unProAvg_0.0_unPiAvg_0.0_ElPro_0.0_InePro_1.0_merge_voxelsIdeal_277e-5m_554e-5m.root");
  //ReconTree->Add("~/TREx/hptrex/src/Data/tobyData/MergedFiles280318/Processed/277_554/TRExRecon_unProAvg_0.0_unPiAvg_0.0_ElPro_1.0_InePro_0.0_merge_voxelsIdeal_277e-5m_554e-5m.root");
  //ReconTree->Add("~/TREx/hptrex/src/Data/tobyData/MergedFiles280318/Processed/277_554/TRExRecon_unProAvg_0.0_unPiAvg_1.0_ElPro_0.0_InePro_0.0_merge_voxelsIdeal_277e-5m_554e-5m.root");
  //ReconTree->Add("~/TREx/hptrex/src/Data/tobyData/MergedFiles280318/Processed/277_554/TRExRecon_unProAvg_1.0_unPiAvg_0.0_ElPro_0.0_InePro_0.0_merge_voxelsIdeal_277e-5m_554e-5m.root");
  //ReconTree->Add("~/TREx/hptrex/src/Data/tobyData/MergedFiles280318/Processed/277_554/TRExRecon_unProAvg_3.0_unPiAvg_3.0_ElPro_0.0_InePro_0.0_merge_voxelsIdeal_277e-5m_554e-5m.root");
  //ReconTree->Add("~/TREx/hptrex/src/Data/tobyData/MergedFiles280318/Processed/277_554/TRExRecon_unProAvg_5.7_unPiAvg_21.7_ElPro_1.0_InePro_1.0_merge_voxelsIdeal_277e-5m_554e-5m.root");
  //ReconTree->Add("~/TREx/hptrex/src/Data/tobyData/MergedFiles280318/Processed/277_554/TRExRecon_unProAvg_6.7_unPiAvg_21.7_ElPro_0.0_InePro_1.0_merge_voxelsIdeal_277e-5m_554e-5m.root");
  //ReconTree->Add("~/TREx/hptrex/src/Data/tobyData/MergedFiles280318/Processed/277_554/TRExRecon_unProAvg_6.7_unPiAvg_21.7_ElPro_1.0_InePro_0.0_merge_voxelsIdeal_277e-5m_554e-5m.root");
  //  ReconTree->Add("~/TREx/hptrex/src/Data/tobyData/MergedFiles280318/Processed/277_554/TRExRecon_unProAvg_7.7_unPiAvg_21.7_ElPro_0.0_InePro_0.0_merge_voxelsIdeal_277e-5m_554e-5m.root");
  ReconTree->Add(inputFile);

  //declare histograms
  TH1D dEdx("dEdx", "dEdx", 50, 0, 0.00002);
  TH1D TrackLength("TrackLength", "TrackLength", 300, 0, 1200);
  TH1D ChargeSum("ChargeSum", "ChargeSum", 200, 0, 0.024);

  TH1D Cleanliness("Cleanliness", "Cleanliness", 100, 0, 1);
  
  TH1D* CleanlinessVsM[21];
  TH1D* CompletenessVsM[21];
  TH2D* TrackEfficiencyScan[21];

  int compbins=10;
  int cleanbins=10;
  double step=0.05;


  for(int i=0; i<21; ++i){    
    char nameClean[20];
    char nameComp[20];
    sprintf(nameClean,"Cleanliness_Mult%d", i);
    sprintf(nameComp, "Completeness_Mult%d", i);
    CleanlinessVsM[i] = new TH1D(nameClean , "", 100, 0, 1);   
    CompletenessVsM[i] = new TH1D(nameComp, "", 100, 0, 1);
    
    char nameScan[50];
    sprintf(nameScan, "CleanlinessVsCompletenessCuts_Mult%d", i);
    TrackEfficiencyScan[i]=new TH2D(nameScan,"", cleanbins, 0.5, 1.00, compbins, 0.5, 1.00);  
  }
  
  TH1D Completeness("Completeness", "Completeness", 100, 0, 1);

  TH1D MeanCleanliness("MeanCleanliness", "MeanCleanliness", 20, 0, 20);
  TH1D MeanCompleteness("MeanCompleteness", "MeanCompleteness", 20, 0, 20);
  
  TH1D TrackEfficiency("TrackEfficiency", "TrackEfficiency", 20, 0, 20);
  TH1D TrackPurity("TrackPurity", "TrackPurity", 20, 0, 20);

  TH1D ProtonEfficiency("ProtonEfficiency", "ProtonEfficiency",20, 0, 20); 
  TH1D ProtonPurity("ProtonPurity", "ProtonPurity", 20,0,20);

  TH2D Multiplicities("Multiplicities", "Multiplicities", 20, 0, 20, 20, 0, 20); 

  TH1D PDGs("PDG", "PDG", 3, 0, 3);
  TH1D GoodPDGs("GoodPDGs", "GoodPDGs", 3, 0, 3);
  TH1D BadPDGs("BadPDGs", "BadPDGs", 3, 0, 3);

  int entries = ReconTree->GetEntries();
  
  trex::WritableEvent * event=0;

  double pion_veto = 3.5e-6;
  double SingleCleanCut = 0.0;
  double SingleCompCut = 0.5;

  ReconTree->SetBranchAddress("event", &event);
   
  //loop over events
  
  int goodProtons[21]={};
  int badProtons[21]={};
  int missedProtons[21]={};
  
  int totalNumberOfPaths[21]={};
  int totalNumberOfTrueTracks[21]={};
  int totalNumberOfGoodTracks[21][10][10]={};

  int totalTracksPassingSingleCut[21]={};
  
  for(int i=0; i<entries; ++i) {
    
    ReconTree->GetEntry(i);

   
    int recoMulti = event->RecoMultiplicity;
    int trueMulti = event->TrueMultiplicity;
    
    //std::cout << "The TRUE MULTIPLICITY was: " << trueMulti << std::endl;
    //std::cout << "The RECO MULTIPLICITY was: " << recoMulti << std::endl;
    
    if (trueMulti >20){
      continue;
    }

    Multiplicities.Fill(trueMulti, recoMulti);
    totalNumberOfTrueTracks[trueMulti]+=trueMulti;


    //Loop over patterns
    std::vector<trex::WritablePattern> pats = event->patterns;
    for(int j=0; j<pats.size(); ++j){
      

      std::vector<double> dEdx_vec = pats[j].dEdx;
      std::vector<double> TrackLength_vec = pats[j].TrackLength;
      std::vector<double> ChargeSum_vec = pats[j].ChargeSum;
      //declare more vectors
      std::vector<double> Cleanliness_vec = pats[j].TrackCleanliness;
      std::vector<double> Completeness_vec = pats[j].TrackCompleteness;
      std::vector<int> PDG_vec = pats[j].PDG;
      std::vector<int> PID_vec = pats[j].PID;
      //std::vector<int> PDG_vec = pats[j].PDG;

      //Loop over paths
      for(int k=0; k<dEdx_vec.size(); ++k){
	
	totalNumberOfPaths[trueMulti]+=1;
	
	dEdx.Fill(dEdx_vec[k]);
	TrackLength.Fill(TrackLength_vec[k]);
	ChargeSum.Fill(ChargeSum_vec[k]);
	Cleanliness.Fill(Cleanliness_vec[k]);
	Completeness.Fill(Completeness_vec[k]);
	

	for(int cleanbin=0; cleanbin<cleanbins; ++cleanbin){
	  
	  for(int compbin=0; compbin<compbins; ++compbin){
	    
	    double compcut = 0.5+compbin*step;
	    double cleancut = 0.5+cleanbin*step;
	    
	    if (Cleanliness_vec[k] > cleancut && Completeness_vec[k] > compcut){
	      totalNumberOfGoodTracks[trueMulti][cleanbin][compbin]+=1;
	    }
	  }
	}
	
	CleanlinessVsM[trueMulti]->Fill(Cleanliness_vec[k]);
	CompletenessVsM[trueMulti]->Fill(Completeness_vec[k]);

	std::cout << "PDG_vec[" << k << "] is = " << PDG_vec[k] << std::endl;
	//apply single cut here
	if(Cleanliness_vec[k] > SingleCleanCut && Completeness_vec[k] > SingleCompCut){

	  //CleanlinessVsM[trueMulti]->Fill(Cleanliness_vec[k]);
	  //CompletenessVsM[trueMulti]->Fill(Completeness_vec[k]); 
	  if(PDG_vec[k]==2212 && PID_vec[k]==-1){
	    goodProtons[trueMulti]+=1;
	  }else if(PDG_vec[k]==211 && PID_vec[k]==-1){
	    badProtons[trueMulti]+=1;
	  }else if(PDG_vec[k]==2212 && PID_vec[k]==1){
	    missedProtons[trueMulti]+=1;
	  }
	  
	  //fill more histograms
	
	  totalTracksPassingSingleCut[trueMulti]+=1;
  
	  //Check PDG here
	  
	  if(PDG_vec[k]==2212){
	    std::cout << "PROTON!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
	  }else if(PDG_vec[k]==211){
	    std::cout << "PION!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
	  }else{std::cout << PDG_vec[k] << std::endl;}

	}
      }      
    }
  } 
  
  
  for(int i=1; i<21; ++i) {
    
    std::cout << "MEAN Cleanliness of multiplicity " << i << " was: " << CleanlinessVsM[i]->GetMean() << std::endl;
    
    if(i==1){
      if (CleanlinessVsM[1]->GetMean()>1.00){
	MeanCleanliness.SetBinContent(1,1);
      }
    }else{
      MeanCleanliness.SetBinContent(i,CleanlinessVsM[i]->GetMean());
    }
    if(i>1){MeanCleanliness.SetBinContent(i,CleanlinessVsM[i]->GetMean());}

    MeanCleanliness.SetBinError(i,CleanlinessVsM[i]->GetRMS());
    MeanCompleteness.SetBinContent(i, CompletenessVsM[i]->GetMean());
    MeanCompleteness.SetBinError(i, CompletenessVsM[i]->GetRMS());

    //Fill the 2D scan for every multiplicity
    for(int cleanbin=1; cleanbin<11; ++cleanbin){
      for(int compbin=1; compbin<11; ++compbin){

	double goodTracks = (double)totalNumberOfGoodTracks[i][cleanbin-1][compbin-1];	
	double trackEfficiency = round(goodTracks/(double)totalNumberOfPaths[i]*100)/100;
	
	TrackEfficiencyScan[i]->SetBinContent(cleanbin, compbin, trackEfficiency);
	
      }
    }
       
    double SingleCutPurity = round((double)totalTracksPassingSingleCut[i]/(double)totalNumberOfPaths[i]*100)/100;
    double SingleCutEfficiency = round((double)totalTracksPassingSingleCut[i]/(double)totalNumberOfTrueTracks[i]*100)/100;
    
    TrackPurity.SetBinContent(i,SingleCutPurity);
    TrackEfficiency.SetBinContent(i,SingleCutEfficiency);    
    
    int numberOfTrueProtons = goodProtons[i]+missedProtons[i];
    int numberOfRecoProtons = goodProtons[i]+badProtons[i];


    std::cout << "Have found " << goodProtons[i] << " good Protons, " << missedProtons[i] << " missed Protons and " << badProtons[i] << " bad Protons." << std::endl; 

    double proEff = (double)goodProtons[i]/(double)numberOfTrueProtons;
    double proPur = (double)goodProtons[i]/(double)numberOfRecoProtons;
    
    ProtonEfficiency.SetBinContent(i,proEff);
    ProtonPurity.SetBinContent(i,proPur);    
    
    std::cout << "Proton Efficiency and Purity was: " << proEff << " and " << proPur << std::endl;
  }

  
  dEdx.SetLineColor(kBlue);
  TrackLength.SetLineColor(kBlue);
  ChargeSum.SetLineColor(kRed);
  
  TFile outf("PID.root", "RECREATE");



  for(int i=1; i<21; ++i) {
    TCanvas *C = new TCanvas("C","",750,750);
    C->cd();
    gStyle->SetOptStat("");
    TrackEfficiencyScan[i]->GetXaxis()->SetTitle("Cleanliness Cut");
    TrackEfficiencyScan[i]->GetYaxis()->SetTitle("Completeness Cut");
    TrackEfficiencyScan[i]->GetYaxis()->SetTitleOffset(1.5);
    TrackEfficiencyScan[i]->GetZaxis()->SetRange(0,1);
    TrackEfficiencyScan[i]->SetOption("COL2 TEXT");
    TrackEfficiencyScan[i]->Write();
    //TrackEfficiencyScan[i]->Draw();
    //char canvas[15];
    //sprintf(canvas, "Scan%d.png", i);
    //C->Print(canvas);
    delete C;
  }

  TH2D AverageTrackEfficiency("AverageTrackEff", "", cleanbins, 0.5, 1.00, compbins, 0.5, 1.00);

  for(int m=3; m<17; ++m){
    AverageTrackEfficiency.Add(TrackEfficiencyScan[m]);    
  }

  AverageTrackEfficiency.Scale(1.00/14.00);  

  //JTH: Have disabled writing the histograms to separate .png files for now as it was
  //causing a crash. Everything should still be in the root file.

  TCanvas *C = new TCanvas("C","",750,750);
  C->cd();
  gStyle->SetOptStat("");
  gStyle->SetPaintTextFormat("4.2f");
  AverageTrackEfficiency.GetXaxis()->SetTitle("Cleanliness Cut");
  AverageTrackEfficiency.GetYaxis()->SetTitle("Completeness Cut");
  AverageTrackEfficiency.GetYaxis()->SetTitleOffset(1.5);
  AverageTrackEfficiency.GetZaxis()->SetRange(0,1);
  AverageTrackEfficiency.SetOption("COL2 TEXT");
  AverageTrackEfficiency.Write();
  AverageTrackEfficiency.Draw();
  //C->Print("AverageTrackEfficiency.png");
  
  TrackPurity.GetYaxis()->SetRangeUser(0,1);
  TrackPurity.SetMarkerStyle(kFullTriangleDown);
  TrackPurity.Draw("P");
  //C->Print("TrackPurityMinimalCut.png");

  TrackEfficiency.GetYaxis()->SetRangeUser(0,1);
  TrackEfficiency.SetMarkerStyle(kFullTriangleDown);
  TrackEfficiency.Draw("P");
  //C->Print("TrackEfficiencyMinimalCut.png");

  ProtonEfficiency.GetYaxis()->SetRangeUser(0,1);
  ProtonEfficiency.SetMarkerStyle(kFullTriangleDown);
  ProtonEfficiency.SetMarkerColor(kRed);
  ProtonEfficiency.Draw("P");
  //C->Print("ProtonEfficiency.png");

  ProtonPurity.GetYaxis()->SetRangeUser(0,1);
  ProtonPurity.SetMarkerStyle(kFullTriangleDown);
  ProtonPurity.SetMarkerColor(kRed); 
  ProtonPurity.Draw("P");
  //C->Print("ProtonPurity.png");

  MeanCleanliness.GetYaxis()->SetRangeUser(0,1);
  MeanCleanliness.SetMarkerStyle(kFullTriangleDown);
  MeanCleanliness.Draw("P");
  //C->Print("MeanCleanliness.png");
  
  MeanCompleteness.GetYaxis()->SetRangeUser(0,1);
  MeanCompleteness.SetMarkerStyle(kFullTriangleDown);
  MeanCompleteness.Draw("P");
  //C->Print("MeanCompleteness.png");

  Multiplicities.Draw("COL2");
  //C->Print("Multiplicities.png");

  dEdx.Write();
  TrackLength.Write();
  ChargeSum.Write();
  Cleanliness.Write();
  Completeness.Write();
  MeanCleanliness.Write();
  MeanCompleteness.Write();
  TrackEfficiency.Write();
  TrackPurity.Write();
  ProtonEfficiency.Write();
  ProtonPurity.Write();
  Multiplicities.Write();
  //write more histograms


}

