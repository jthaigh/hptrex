#ifndef TSIMLOADERCCD_HXX
#define TSIMLOADERCCD_HXX


//C++
#include <vector>

// ROOT
#include <TVector3.h>
#include "TFile.h"
#include "TTree.h"
#include <THnSparse.h>
#include "TH3D.h"
#include "TImage.h"
#include "TCanvas.h"

//TREx
#include "TTPCHitPad.hxx"
#include "TTrueHit.hxx"
#include "TTRExPattern.hxx"
#include "TTrueTrack.hxx"

//Sim data
#include "GasTPCDataLib.hxx"

namespace trex{

  class TSimLoader{

  public:

    TSimLoader(std::string inputFile);
    
    void LoadEvent(unsigned int i);
    
    void DrawDetector(){
      TImage *img = TImage::Create();
      TCanvas * c = new TCanvas;
      Detector->SetMarkerColor(kRed);
      Detector->GetXaxis()->SetTitle("X");
      Detector->GetYaxis()->SetTitle("Y");
      Detector->GetZaxis()->SetTitle("Z");
      c->cd();
      Detector->Draw("");
      img->FromPad(c);
      img->WriteImage("Detector.png");
      Detector->Write("Detector");
      delete img;
      delete c;
    }

    unsigned int GetNEvents();


    inline std::vector<trex::TTPCHitPad*>& GetHits(){return fHits;}
    
    inline std::vector<TTrueHit*>& GetTrueHits(){return fTrueHits;}
    
    unsigned int GetNVoxels();
    
    int GetTrueMultiplicity(){return fBP_NParticles;}

    
    struct voxel {
      double x_pos;
      double y_pos;
      double z_pos;
      double time;
      
      double Edep;
      
      void printVoxel() {
	if(Edep<10){
	  std::cout  << "\n" <<  x_pos << " : " << y_pos << " :  " << z_pos << "\n" << "     ---------\n" << "     |       |\n" << "     |   " << Edep << "   |\n" << "     |       |\n" << "     ---------" << std::endl; }
	else{std::cout << "\n" <<  x_pos << " : " << y_pos << " :  " << z_pos << "\n" << "     ---------\n" << "     |       |\n" << "     |  " << Edep << "   |\n" << "     |       |\n" << "     ---------" << std::endl;}
      };     
    };
    
    inline std::vector<voxel*>& GetVoxels(){return fVoxels;}
    

    TFile* GetFile() {
      return fFile;
    }
     
  private:
    
    TTree* fVoxelsTree;
    TTree* fTruthTree;
    TH3D * Detector;
    TFile* fFile;

    std::vector<trex::TTPCHitPad*> fHits;
    std::vector<voxel*> fVoxels;
    std::vector<TTrueHit*> fTrueHits;

    //int fimageNumber = -1;
    std::vector<trex::TTrueTrack*> fTrueTracks;

    //Variables to point to branches in input trees
    THnSparseF* fBP_voxels;
    THnSparseS* fBP_truthmatch_voxels;

    Double_t fBP_Momentum[200];
    Double_t fBP_Xi[200];
    Double_t fBP_Yi[200];
    Double_t fBP_Zi[200];
    Double_t fBP_Xf[200];
    Double_t fBP_Yf[200];
    Double_t fBP_Zf[200];
    Int_t fBP_pdg[200];
    Int_t fBP_TrackID[200];
    Int_t fBP_ParentID[200];
    Int_t fBP_NParticles;

    
  };
}

#endif
