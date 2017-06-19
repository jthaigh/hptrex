
#include <THnSparse.h>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TVector3.h"

#include "../LinkDef.hh"

std::string base_name(std::string const & path);

int  Voxelize(const char * inputfile, int x_voxelDim, int y_voxelDim, Int_t events, string mode) //Define voxelDim in tens of microns i.e. 1e-5m
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
  
  int x_Bins = (int)(1200/x_res);
  int y_Bins = (int)(1200/y_res);
  int z_Bins = 1;
  
  double x_range = x_Bins*x_res;
  double y_range = y_Bins*y_res;

  //Int_t x_res = 1200/x_res;
  //Int_t y_res = 1200/y_res;
  
  //Double_t maxs[dim] = { x_range/2, y_range/2, 1.};
  //Double_t mins[dim] = { -x_range/2, -y_range/2, 0.};
  //Int_t bins[dim] = {x_Bins, y_Bins, z_Bins};
  
  //use original resolution of 2.34 mm
  Double_t maxs[dim] = { 600.21, 600.21, 1.};
  Double_t mins[dim] = { -600.21, -600.21, 0.};
  Int_t bins[dim] = {513, 513, 1};

  
  //setting up new Voxels Tree as TREx input
  TTree* VoxelsTree = new TTree("VoxelsTree", "VoxelsTree");
  THnSparseF* voxels = new THnSparseF("Voxels","", dim, bins, mins, maxs);
  VoxelsTree->Branch("voxels", "THnSparseF", &voxels);  

  //DEFINE TRACK LEVEL TRUTH INFORMATION VARIABLES HERE!
  Int_t ImageNumber;
  VoxelsTree->Branch("ImageNumber", &ImageNumber, "ImageNumber/I");
  Double_t Momentum;
  VoxelsTree->Branch("Momentum", &Momentum, "Momentum/D");
  TVector3 TrueXi; 
  VoxelsTree->Branch("Xi", "TVector3", &TrueXi);
  TVector3 TrueXf;
  VoxelsTree->Branch("Xf", "TVector3", &TrueXf);
  Int_t PDG;
  VoxelsTree->Branch("PDG", &PDG, "PDG.I");
  Int_t TrackID;
  VoxelsTree->Branch("TrackID", &TrackID, "TrackI/I");
  Int_t ParentID;
  VoxelsTree->Branch("ParentID", &ParentID, "ParentID/I");
  Int_t ProOrPi;
  VoxelsTree->Branch("ProOrPi", &ProOrPi, "ProOrPi/I");


  //LOOP FOR READING OUT TTREE RATHER THAN IMAGES
  
  TTree* MergeTree = (TTree*)f.Get("MergeTree");
  TH2D* SingleHist = new TH2D();   

  //Set all Branch Addresses
  if(mode == "Ideal"){
    MergeTree->SetBranchAddress("hIdeal", &SingleHist);}
  else if(mode == "Real"){
    MergeTree->SetBranchAddress("hReal", &SingleHist);}
  else{std::cout << "PLEASE SET A VALID INPUT MODE - EITHER 'Ideal' OR 'Real'" << std::endl;
    return 1;}


  MergeTree->SetBranchAddress("ImageNumber", &ImageNumber);
  MergeTree->SetBranchAddress("P", &Momentum);
  Double_t Xi;
  MergeTree->SetBranchAddress("Xi", &Xi);
  Double_t Yi;
  MergeTree->SetBranchAddress("Yi", &Yi);
  Double_t Zi;
  MergeTree->SetBranchAddress("Zi", &Zi);
  Double_t Xf;
  MergeTree->SetBranchAddress("Xf", &Xf);
  Double_t Yf;
  MergeTree->SetBranchAddress("Yf", &Yf);
  Double_t Zf;
  MergeTree->SetBranchAddress("Zf", &Zf);
  MergeTree->SetBranchAddress("pdg", &PDG);
  MergeTree->SetBranchAddress("TrackID", &TrackID);
  MergeTree->SetBranchAddress("ParentID", &ParentID);
  MergeTree->SetBranchAddress("ProOrPi", &ProOrPi);
 

  
  int entries = MergeTree->GetEntries();

  //Set how many events we want to process
  int nEvents;  

  if(events < entries){
    nEvents=events;
  }else{nEvents=entries;}

  std::cout << "This Tree contains " << entries << " Events." << std::endl;


  //EventLoop
  for(int i=0; i<nEvents; ++i){
    
    std::cout << "Voxelising Event # " << i << std::endl;

    MergeTree->GetEntry(i);
    
    voxels->Reset();

       
    //Reading out Hits
    for (int x=0; x<501; ++x){
      for(int y=0; y<501; ++y){
	
	Double_t Edep = SingleHist->GetBinContent(x,y);
	Double_t xpos = SingleHist->GetXaxis()->GetBinCenter(x);
	Double_t ypos = SingleHist->GetYaxis()->GetBinCenter(y);
	Double_t zpos = 1;

	Double_t position[3] = {xpos, ypos, zpos};
        
	
	if(Edep!=0){
	  voxels->Fill(position, Edep);	
	}
      }
    } 
     
    TrueXi.SetXYZ(Xi, Yi, Zi);
    TrueXf.SetXYZ(Xf, Yf, Zf);
 
    VoxelsTree->Fill();

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
  
  
