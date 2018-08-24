
#include <THnSparse.h>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TVector3.h"
#include "TROOT.h"
#include "../LinkDef.hh"
#include "TObjArray.h"
#include "THnSparse.h"

std::string base_name(std::string const & path);

int  Voxelize(const char * inputfile, int x_voxelDim, int y_voxelDim, Int_t events, string mode, double threshold=0.) //Define voxelDim in tens of microns i.e. 1e-5m. events is number of events to be processed. mode is one of "Ideal" or "Real". 

{

  //////////////////////////////////////////////
  //Twiddling with filenames and opening files//
  //////////////////////////////////////////////

  char inname[100];
  char outname[500];
  char N[100];
  
  std::cout << "The input file is " << inputfile << "\n" << std::endl;
  
  TFile f(inputfile);
  if (f.GetNkeys()<1) {
    printf("Error: Cannot open file\n");
    return 1;
  }  
    
  //Create a sensibly named output file
  string input_path = std::string(inputfile);
  string basefile = base_name(input_path)+"_voxels%s_%de-5m_%de-5m_%g.root";
  const char * basefile_char = basefile.c_str();
  const char * input_mode = mode.c_str();
  sprintf(outname, basefile_char, input_mode, x_voxelDim,y_voxelDim,threshold);
  
  std:: cout << "New file will be called: " << outname << std::endl;
  TFile f1(outname,"RECREATE"); 


  
  //////////////////////////////////////////////
  //Define voxel binning - currently hardcoded//
  //////////////////////////////////////////////

  const Int_t dim = 3;
  
  Double_t maxs[dim] = { 710, 710, 1.};
  Double_t mins[dim] = { -710, -710, 0.};

  Int_t bins[dim] = {1529, 1529, 1};

  /*
  //JTH: This isn't currently used since binning was changed to hardcoding despite
  //parameters being present in macro sig. Misleading, should be addressed ASAP.
  double x_res = x_voxelDim/100; //to get resolution in mm from 10s of microns
  double y_res = y_voxelDim/100;
  
  int x_Bins = (int)(1420/x_res);
  int y_Bins = (int)(1420/y_res);
  int z_Bins = 1;
  
  double x_range = x_Bins*x_res;
  double y_range = y_Bins*y_res;
  */

  //use this for user defined resolution from input
  //Double_t maxs[dim] = { x_range/2, y_range/2, 1.};
  //Double_t mins[dim] = { -x_range/2, -y_range/2, 0.};
  //Int_t bins[dim] = {x_Bins, y_Bins, z_Bins};
  
  //use original resolution of 2.34 mm
  //Double_t maxs[dim] = { 600.21, 600.21, 1.};
  //Double_t mins[dim] = { -600.21, -600.21, 0.};



  /////////////////////////////////////////////////////////////
  //Creation of output VoxelsTree to store converted image   //
  //JTH - Added an extra histogram truthmatch_voxels to store//
  //the information on which track caused the voxel deposit  //
  /////////////////////////////////////////////////////////////
  
  //setting up new Voxels Tree as TREx input
  TTree* VoxelsTree = new TTree("VoxelsTree", "VoxelsTree");

  THnSparseF* voxels = new THnSparseF("Voxels","", dim, bins, mins, maxs);
  VoxelsTree->Branch("voxels", "THnSparseF", &voxels);  

  THnSparseS* truthmatch_voxels = new THnSparseS("truthmatch_Voxels","", dim, bins, mins, maxs);
  VoxelsTree->Branch("truthmatch_voxels", "THnSparseS", &truthmatch_voxels);  

  Int_t EventNumber;
  VoxelsTree->Branch("EventNumber", &EventNumber, "EventNumber/I");



  ////////////////////////////////////////////////////////////////
  //Creation of output TruthTree to store track-level truth info//
  //from input file.                                            //
  ////////////////////////////////////////////////////////////////

  //JTH TODO: Much better to have 1 entry in the TruthTree per event rather
  //than one per track, and store particle-level quantities in arrays/vectors.
  //However this requires changing main TREx TSimLoader code so will defer for now.

  TTree* TruthTree = new TTree("TruthTree","TruthTree");

  /*  Double_t MOMENTUM;
  TruthTree->Branch("Momentum", &MOMENTUM, "Momentum/D");

  TVector3 TrueXi; 
  TruthTree->Branch("Xi", "TVector3", &TrueXi);

  TVector3 TrueXf;
  TruthTree->Branch("Xf", "TVector3", &TrueXf);

  Int_t TRACKID;
  TruthTree->Branch("TrackID", &TRACKID, "TrackI/I");

  Int_t PARENTID;
  TruthTree->Branch("ParentID", &PARENTID, "ParentID/I");

  Int_t PDG;
  TruthTree->Branch("pdg", &PDG, "pdg/I");

  Int_t NPARTICLES;
  TruthTree->Branch("NParticles", &NPARTICLES, "NParticles/I");
  */

  ////////////////////////////////////////////////////////////////////////
  //Get event-level image/truth tree from input file and assign branches//
  //to variables.                                                       //
  ////////////////////////////////////////////////////////////////////////

  TTree* MergeTree = (TTree*)f.Get("MergeTree");
 

  //Histograms of composite event-level images (including noise in Real case).
  TH2D* MergeHist = 0;// = new TH2D();

  //Histograms of individual tracks.
  TObjArray* SingleHist = 0;

  //JTH: Not needed?
  //  SingleHist = new TObjArray(); //Create the TH2D array 
  //SingleHist->SetOwner(kTRUE);

  if(mode == "Ideal"){
    MergeTree->SetBranchAddress("EventImageIdeal", &MergeHist);
    MergeTree->SetBranchAddress("ImageArrayIdeal", &SingleHist);
  }
  else if(mode == "Real"){
    MergeTree->SetBranchAddress("EventImageReal", &MergeHist);
    MergeTree->SetBranchAddress("ImageArrayReal", &SingleHist);
  }
  else{std::cout << "PLEASE SET A VALID INPUT MODE - EITHER 'Ideal' OR 'Real'" << std::endl;
    return 1;}

  MergeTree->SetBranchAddress("EventNumber", &EventNumber);
  Double_t Momentum[200];
  MergeTree->SetBranchAddress("P", Momentum);
  Double_t Xi[200];
  MergeTree->SetBranchAddress("Xi", Xi);
  Double_t Yi[200];
  MergeTree->SetBranchAddress("Yi", Yi);
  Double_t Zi[200];
  MergeTree->SetBranchAddress("Zi", Zi);
  Double_t Xf[200];
  MergeTree->SetBranchAddress("Xf", Xf);
  Double_t Yf[200];
  MergeTree->SetBranchAddress("Yf", Yf);
  Double_t Zf[200];
  MergeTree->SetBranchAddress("Zf", Zf);
  Int_t pdg[200];
  MergeTree->SetBranchAddress("pdg", pdg);
  Int_t TrackID[200];
  MergeTree->SetBranchAddress("TrackID", TrackID);
  Int_t ParentID[200];
  MergeTree->SetBranchAddress("ParentID", ParentID);
  Int_t NParticles;
  MergeTree->SetBranchAddress("NParticles",&NParticles);

  //Make branches in truth tree to pass through truth info
  //from input file.
  TruthTree->Branch("NParticles", &NParticles, "NParticles/I");
  TruthTree->Branch("Momentum", &Momentum, "Momentum[NParticles]/D");
  TruthTree->Branch("Xi", Xi, "Xi[NParticles]/D");
  TruthTree->Branch("Yi", Yi, "Yi[NParticles]/D");
  TruthTree->Branch("Zi", Zi, "Zi[NParticles]/D");
  TruthTree->Branch("Xf", Xf, "Xf[NParticles]/D");
  TruthTree->Branch("Yf", Yf, "Yf[NParticles]/D");
  TruthTree->Branch("Zf", Zf, "Zf[NParticles]/D");
  TruthTree->Branch("TrackID", TrackID, "TrackID[NParticles]/I");
  TruthTree->Branch("ParentID", ParentID, "ParentID[NParticles]/I");
  TruthTree->Branch("pdg", pdg, "pdg[NParticles]/I");
  
  //####################################//
  //Main event loop                     //
  //####################################//

  int entries = MergeTree->GetEntries();
  std::cout << "This Tree contains " << entries << " Spills." << std::endl;
  for(int i=0; i<entries&&i<events; ++i){

    ///////////////////////////////////
    //Book-keeping and console output//
    ///////////////////////////////////

    MergeTree->GetEntry(i);

    std::cout << "Voxelising Spill # " << i << std::endl;
    std::cout << "NParticles = " << NParticles << std::endl;

    voxels->Reset();
    truthmatch_voxels->Reset();

    TruthTree->Fill();

    ///////////////////////////////////////////////
    //Pass through per-particle truth information//
    //and fill the truth-matching histogram      //
    ///////////////////////////////////////////////
    //JTH: If multiple input bins correspond to a single voxel then this book-keeping logic
    //will not function properly since I am mapping the bins in the input THnSparse
    //straight onto the output histogram. Would require a modification to the code.

    //Temporary histogram to store the value of the biggest single-particle
    //contribution to a voxel.
    TH2D truthmatch_biggest_contributor("","", bins[0], mins[0], maxs[0],bins[1], mins[1], maxs[1]);
    truthmatch_biggest_contributor.SetDirectory(0);
    
    for(int j=0; j<NParticles; ++j){
      std::cout << "Particle number = " << j << " of " << NParticles << std::endl;
      
      //Per-particle truth variables
      /*      TrueXi.SetXYZ(Xi[j], Yi[j], Zi[j]);    
      TrueXf.SetXYZ(Xf[j], Yf[j], Zf[j]);
      PDG = pdg[j];
      TRACKID = TrackID[j];
      PARENTID = ParentID[j];
      MOMENTUM = Momentum[j];
      NPARTICLES = NParticles;
      */
 
      Int_t coord[3];
      coord[2]=0;
      THnSparse* sparseHist=(THnSparse*) SingleHist->At(j);
      Long64_t nBins=sparseHist->GetNbins();
      int countFills=0;
      for(int iBin=0;iBin<nBins;++iBin){
	double Edep = sparseHist->GetBinContent(iBin,coord);
	if(truthmatch_biggest_contributor.GetBinContent(coord[0],coord[1])<Edep){
	  truthmatch_biggest_contributor.SetBinContent(coord[0],coord[1],Edep);
	  truthmatch_voxels->SetBinContent(coord,j+1);
	  ++countFills;
	}
      }
      std::cout<<"Particle had "<<nBins<<" non-zero bins in image of which we used "<<countFills<<std::endl;
    }

    //Reading out Hits
    for (int x=0; x<1529; ++x){
      for(int y=0; y<1529; ++y){
	//	Double_t xpos = MergeHist->GetXaxis()->GetBinCenter(x);
	//Double_t ypos = MergeHist->GetYaxis()->GetBinCenter(y);
	//Double_t zpos = 0.5;
	Double_t Edep = MergeHist->GetBinContent(x,y);
	//Double_t position[3] = {xpos, ypos, zpos};
	
	if(Edep>threshold){
	  Int_t position[3] = {x,y,0};
	  voxels->SetBinContent(position, Edep);	
	}
      }
    }
    
    VoxelsTree->Fill();
  }
  
  
  std::cout << "All Events written to Voxel Tree" << std::endl;
  
  TruthTree->Print();
  
  std::cout << "TruthTree printed successfully " << std::endl;

  VoxelsTree->Print();
  
  std::cout << "VoxelTree printed succesfully" << std::endl;
  
  std:: cout << "OUTPUT CAN BE FOUND HERE: " << outname << std::endl;
  
  // Write Output to File
  f1.Write();
  f1.Close();  
  return 0;
    
}

  
  std::string base_name(std::string const & path)
  {
    string base =  path.substr(path.find_last_of("/\\") + 1);
    
    return base.substr(0,base.length()-5);
    
    //return path.substr(0, path.length()-5);
  }
  
  
