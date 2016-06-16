#ifndef TSIMLOADER_HXX
#define TSIMLOADER_HXX

//C++
#include <vector>

// ROOT
#include <TVector3.h>
#include "TFile.h"
#include "TTree.h"
#include <THnSparse.h>

//TREx
#include "TTPCHitPad.hxx"

//Sim data
#include "GasTPCDataLib.hxx"

namespace trex{

  class TSimLoader{

  public:

    TSimLoader(std::string inputFile);

    void LoadEvent(unsigned int i);

    unsigned int GetNEvents();


    inline std::vector<trex::TTPCHitPad*>& GetHits(){return fHits;}

    
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


  private:

    TFile* fFile;

    TTree* fTree;
    TTree* fVoxelsTree;

    std::vector<trex::TTPCHitPad*> fHits;
    std::vector<voxel*> fVoxels;

    SimulData* fSimulDataBranch;
    THnSparseF* fVoxelBranch;

  };
}

#endif
