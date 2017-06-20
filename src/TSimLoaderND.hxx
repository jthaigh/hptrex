#ifndef TSIMLOADERND_HXX
#define TSIMLOADERND_HXX


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

//Sim data
#include "GasTPCDataLib.hxx"
#include "evtrec/EventRecord.h"
#include "evtrec/MCHit.h"
#include "evtrec/Vector3D.h"
#include "evtrec/Vector4D.h"
#include "evtrec/Units.h"
//#include "evtrec/RootFileWriter.h"
//#include "evtrec/RootFileReader.h"
#include "evtrec/RecoParticle.h"
#include "evtrec/MCTrack.h"
#include "evtrec/MCParticle.h"
//#include "evtrec/MCGenInfo.h"
#include "evtrec/LinkDef.h" 

namespace trex{

  class TSimLoaderND{

  public:

    TSimLoaderND(std::string inputFile);
    
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
    inline std::vector<TTrueTrack*>& GetTrueTracks(){return fTrueTracks;}
    
    unsigned int GetNVoxels();
    
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

    //TTree* GetReconTree(){
    //return fReconTree;
    //}

    TFile* fFile;

     
  private:
    
    //TFile* fFile;
    TFile* fReconOutFile;
    TTree* fTree;
    //TTree* fVoxelsTree;
    //TTree* fReconTree;
    TH3D * Detector;
    
    std::vector<trex::TTPCHitPad*> fHits;
    std::vector<voxel*> fVoxels;
    std::vector<TTrueHit*> fTrueHits;
    std::vector<trex::TTrueTrack*> fTrueTracks;
    
    gastpc::EventRecord* fEventRecordBranch;
    //GeantTrackingTruth* fGeantBranch;
    //THnSparseF* fVoxelBranch;
    
  };
}

#endif
