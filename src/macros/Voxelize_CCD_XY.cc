
#include <THnSparse.h>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TVector3.h"
#include "TROOT.h"
#include "../LinkDef.hh"
#include "TObjArray.h"


std::string base_name(std::string const & path);

<<<<<<< HEAD
int  Voxelize(const char * inputfile, int x_voxelDim, int y_voxelDim, Int_t events, string mode)  //Define voxelDim in tens of microns i.e. 1e-5m
=======


int  Voxelize(const char * inputfile, int x_voxelDim, int y_voxelDim, Int_t events, string mode) //Define voxelDim in tens of microns i.e. 1e-5m. events is number of events to be processed. mode is one of "Ideal" or "Real". 
>>>>>>> 6f00b7874a1b0801b88832bb262700781aae568e
{

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
  string basefile = base_name(input_path)+"_voxels%s_%de-5m_%de-5m.root";
  const char * basefile_char = basefile.c_str();
  const char * input_mode = mode.c_str();
  sprintf(outname, basefile_char, input_mode, x_voxelDim,y_voxelDim);
  
  std:: cout << "New file will be called: " << outname << std::endl;
  TFile f1(outname,"RECREATE"); 
  
  //Setup voxel histogram dimensions with resolution defined by user
  const Int_t dim = 3;
  
  double x_res = x_voxelDim/100; //to get resolution in mm from 10s of microns
  double y_res = y_voxelDim/100;
  
  int x_Bins = (int)(1420/x_res);
  int y_Bins = (int)(1420/y_res);
  int z_Bins = 1;
  
  double x_range = x_Bins*x_res;
  double y_range = y_Bins*y_res;
  
  //use this for user defined resolution from input
  //Double_t maxs[dim] = { x_range/2, y_range/2, 1.};
  //Double_t mins[dim] = { -x_range/2, -y_range/2, 0.};
  //Int_t bins[dim] = {x_Bins, y_Bins, z_Bins};
  
<<<<<<< HEAD
  //use original resolution of 2.34 mm
  //Double_t maxs[dim] = { 600.21, 600.21, 1.};
  //Double_t mins[dim] = { -600.21, -600.21, 0.};
  
  Double_t maxs[dim] = { 710, 710, 1.};
  Double_t mins[dim] = { -710, -710, 0.};
=======
  //use this for original resolution of 2.34 mm
  Double_t maxs[dim] = { 600.21, 600.21, 1.};
  Double_t mins[dim] = { -600.21, -600.21, 0.};
>>>>>>> 6f00b7874a1b0801b88832bb262700781aae568e
  Int_t bins[dim] = {513, 513, 1};
  
  //setting up new Voxels Tree as TREx input
  TTree* VoxelsTree = new TTree("VoxelsTree", "VoxelsTree");
  THnSparseF* voxels = new THnSparseF("Voxels","", dim, bins, mins, maxs);
  VoxelsTree->Branch("voxels", "THnSparseF", &voxels);  

  //DEFINE TRACK LEVEL TRUTH INFORMATION VARIABLES HERE!
  Int_t EventNumber;
  VoxelsTree->Branch("EventNumber", &EventNumber, "EventNumber/I");
  Double_t MOMENTUM;
  VoxelsTree->Branch("Momentum", &MOMENTUM, "Momentum/D");
  TVector3 TrueXi; 
  VoxelsTree->Branch("Xi", "TVector3", &TrueXi);
  TVector3 TrueXf;
  VoxelsTree->Branch("Xf", "TVector3", &TrueXf);
//  Int_t PDG;
//  VoxelsTree->Branch("PDG", &PDG, "PDG.I");
  Int_t TRACKID;
  VoxelsTree->Branch("TrackID", &TRACKID, "TrackI/I");
  Int_t PARENTID;
  VoxelsTree->Branch("ParentID", &PARENTID, "ParentID/I");
  Int_t PDG;
  VoxelsTree->Branch("pdg", &PDG, "pdg/I");
  									// Possible add NParticles to this tree?


  //LOOP FOR READING OUT TTREE RATHER THAN IMAGES
  
  TTree* MergeTree = (TTree*)f.Get("MergeTree");
 // TH2D* SingleHist = new TH2D();   					//Object arrays or TH2D? 
  TObjArray* SingleHist;
  SingleHist = new TObjArray(); //Create the TH2D array                                                                                                      
  SingleHist->SetOwner(kTRUE);


  //Set all Branch Addresses
  if(mode == "Ideal"){
    MergeTree->SetBranchAddress("ImageArrayIdeal", &SingleHist);}
  else if(mode == "Real"){
    MergeTree->SetBranchAddress("ImageArrayReal", &SingleHist);}
  else{std::cout << "PLEASE SET A VALID INPUT MODE - EITHER 'Ideal' OR 'Real'" << std::endl;
    return 1;}

  MergeTree->SetBranchAddress("EventNumber", &EventNumber);
  Double_t Momentum[200];
  MergeTree->SetBranchAddress("P", &Momentum);
  Double_t Xi[200];
  MergeTree->SetBranchAddress("Xi", &Xi);
  Double_t Yi[200];
  MergeTree->SetBranchAddress("Yi", &Yi);
  Double_t Zi[200];
  MergeTree->SetBranchAddress("Zi", &Zi);
  Double_t Xf[200];
  MergeTree->SetBranchAddress("Xf", &Xf);
  Double_t Yf[200];
  MergeTree->SetBranchAddress("Yf", &Yf);
  Double_t Zf[200];
  MergeTree->SetBranchAddress("Zf", &Zf);
  Int_t pdg[200];
  MergeTree->SetBranchAddress("pdg", &pdg);
  Int_t TrackID[200];
  MergeTree->SetBranchAddress("TrackID", &TrackID);
  Int_t ParentID[200];
  MergeTree->SetBranchAddress("ParentID", &ParentID);
  Int_t NParticles;
  MergeTree->SetBranchAddress("NParticles",&NParticles);


  int entries = MergeTree->GetEntries();

  //Set how many events we want to process
  int nEvents;  

  if(events==0){
    nEvents=entries;
  }else if(events < entries){
    nEvents=events;
  }else{nEvents=entries;}

  std::cout << "This Tree contains " << entries << " Spills." << std::endl;


  TH2D* testHist = new TH2D();


  //EventLoop
  for(int i=0; i<nEvents; ++i){

    //voxels->Reset();

    std::cout << "Voxelising Spill # " << i << std::endl;
    
    MergeTree->Print();

    MergeTree->GetEntry(i);

    std::cout << "NParticles = " << NParticles << std::endl;
    
    for(int j=0; j<NParticles; ++j){

    	std::cout << "Particle number = " << j << " of " << NParticles << std::endl;

    	voxels->Reset();

    	//Reading out Hits
    	for (int x=0; x<513; ++x){
      		for(int y=0; y<513; ++y){
	
			testHist = (TH2D*)SingleHist->At(j); 					//->GetBinContent(x,y);
			Double_t xpos = testHist->GetXaxis()->GetBinCenter(x);
			Double_t ypos = testHist->GetYaxis()->GetBinCenter(y);
			Double_t zpos = 1;
			Double_t Edep = testHist->GetBinContent(x,y);
			Double_t position[3] = {xpos, ypos, zpos};
	
			if(Edep!=0){
	  		voxels->Fill(position, Edep);	
			}
      		}
    	} 
        
	//Fill Voxel Tree
    	TrueXi.SetXYZ(Xi[j], Yi[j], Zi[j]);    
    	TrueXf.SetXYZ(Xf[j], Yf[j], Zf[j]);
    	PDG = pdg[j];
	TRACKID = TrackID[j];
	PARENTID = ParentID[j];
	MOMENTUM = Momentum[j];
	VoxelsTree->Fill();
    
    }
    //VoxelsTree->Fill();
  }
  

  /*
  //Reading out 2D images - obsolete?
  
  for(Int_t h=0; h<events; ++h){
    
    voxels->Reset();
    
    //Use user-defined input mode - possible values "Ideal" or "Real" 
    sprintf(inname,"hReadOut%s_img%.3d", input_mode, h);

    std::string files;
    
    //std::cout << "Histogram being read is: " << inname << std::endl;
    
    sprintf(N, "histo%.3d", h);
    
    TH2D *N=(TH2D*)f.Get(inname);
    
    Double_t minY = N->GetBin(0,0);
    Double_t minX = N->GetBin(0,0);
    Double_t maxX = N->GetBin(512,512);
    Double_t maxY = N->GetBin(512,512);    
    
    for(int i=1; i<501; ++i){
      
      for(int j=1; j<501; ++j){
	
	Double_t Edep = N->GetBinContent(i,j);
	
	if(Edep!=0){
	  
	  Double_t xpos = i*2+mins[0];
	  Double_t ypos = j*2+mins[1];
	  Double_t zpos = 1;
	  
	  Double_t position[3] = {xpos, ypos, zpos};
	  
	  //std::cout << "Filling hist with voxel at : " << position[0] << ", " << position[1] << ", " << position[2] << std::endl;
	  voxels->Fill(position, Edep);
	}
      }
    }
  }*/

  std::cout << "All Events written to Voxel Tree" << std::endl;
  
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
    //string base =  path.substr(path.find_last_of("/\\") + 1);
    
    //return base.substr(0,base.length()-5);
    
    return path.substr(0, path.length()-5);
  }
  
  
