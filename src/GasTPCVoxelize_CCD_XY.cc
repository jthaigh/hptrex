
#include <THnSparse.h>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"

#include "GasTPCDataLib.hxx"
#include "LinkDef.hh"

int  GasTPCVoxelize(const char * inputfile, Int_t x_voxelDim, Int_t y_voxelDim, Int_t NFiles) //Define voxelDim in mm
{
  
  char inname[100];
  char outname[100];
  char N[100];
  
  std::cout << "The input file is " << inputfile << "\n" << std::endl;

//work in progress
  //TFile f(inputfile.c_str(),"UPDATE");  //input file from screen
  TFile f(inputfile);
  if (f.GetNkeys()<1) {
    printf("Error: Cannot open file\n");
    return 1;
  }  

 /* TTree *TrackingTree = (TTree*) f.Get("Tracking");  			 //old
  TrackingTree->SetBranchStatus("event",0);
  SimulData *SimulData = 0;
  TrackingTree->SetBranchAddress("SimulData", &SimulData);*/

////////////////////////////////////////////////////////////////////////////////////////////
 const Int_t dim = 3;

 sprintf(outname, "voxelresults_%dmm_%dmm.root", x_voxelDim,y_voxelDim);
 std:: cout << "New file will be called: " << outname << std::endl;
 TFile f1(outname,"RECREATE"); 

 // Int_t res = 1024/voxelDim;
 //P.D.
 Int_t x_res = 1000/x_voxelDim;
 Int_t y_res = 1000/y_voxelDim;
 Int_t res = 1;

 Double_t maxs[dim] = { 500., 500., 1};
 Double_t mins[dim] = { -500., -500., 0};
 Int_t bins[dim] = {x_res, y_res, res};


 THnSparseF* voxels = new THnSparseF("Voxels","", dim, bins, mins, maxs);


 TTree* VoxelsTree = new TTree("VoxelsTree", "VoxelsTree");
 VoxelsTree->Branch("voxels", "THnSparseF", &voxels);


 for(Int_t h=0; h<NFiles; ++h){

  voxels->Reset();

  string mode = "Ideal";
  const char * input_mode = mode.c_str();
  sprintf(inname,"hReadOut%s_img%.3d", input_mode, h);
  //sprintf(inname,"hReadOutIdeal_img%.3d", h);
  //sprintf(inname,"hReadOutReal_img%.3d", h);
  std::string files;

  std::cout << "Histogram being read is: " << inname << std::endl;

  sprintf(N, "histo%.3d", h);

  TH2D *N=(TH2D*)f.Get(inname);
  
  Double_t minY = N->GetBin(0,0);
  Double_t minX = N->GetBin(0,0);
  Double_t maxX = N->GetBin(512,512);
  Double_t maxY = N->GetBin(512,512);

  std::cout << "Min X = " << minX << std::endl;
  std::cout << "Max X = " << maxX << std::endl;
  std::cout << "Bins = " << maxX - minX << std::endl;




/*  Double_t maxs[dim] = { maxX, maxY, 512.};    //need to think about bins
  Double_t mins[dim] = { minX, minY, 512.};
  //Int_t bins[dim] = {res, res, res};
  THnSparseF* voxels = new THnSparseF("Voxels","", dim, bins, mins, maxs); */


  for(int i=1; i<501; ++i){

    for(int j=1; j<501; ++j){

      Double_t Edep = N->GetBinContent(i,j);
      
      if(Edep!=0){
	
	Double_t xpos = i*2+mins[0];
	Double_t ypos = j*2+mins[1];
	Double_t zpos = 1;
	
	Double_t position[3] = {xpos, ypos, zpos};
	
	std::cout << "Filling hist with voxel at : " << position[0] << ", " << position[1] << ", " << position[2] << std::endl;
	voxels->Fill(position, Edep);
      }
    }
  }
  


/*
  for(Double_t b1=1.; b1<513.; ++b1){
	
	for(Double_t b2=1.; b2<513.; ++b2){

		Double_t bincont = N->GetBinContent(b1,b2);  //bin content is the E deposition?

		if (bincont != 0){
			Double_t xyz[3] = {((b1*2.)-500.-(12.*512./b1)),((b2*2.)-500.-(12.*512./b2)),1.};
			std::cout << "Bin Content of " << b1 << "," << b2 << " is " << bincont << std::endl;
			voxels->Fill(xyz,bincont*1000);
		}

	}

  }
*/
 
 VoxelsTree->Fill();

 }


////////////////////////////////////////////////////////////////////////////////////////////
 
//  std::cout << "\nLoaded trees succesfully\n" << std::endl;

// Create File
/*  
  char outname[100];
  sprintf(outname, "Results/g4_nu_190415/Voxel/um/voxelresults_%dum.root", voxelDim);
  std:: cout << "New file will be called: " << outname << std::endl;
  TFile f1(outname,"RECREATE"); 

 
//Create Tree
  TTree* VoxelsTree = new TTree("VoxelsTree", "VoxelsTree");

  std::cout << "New trees created\n" << std::endl;

  const Int_t dim = 3;

//  Int_t bins[dim] = {700, 700, 700};

  

  Int_t res = 7000000/voxelDim; // Convert input (e.g. mm, um, cm, etc.) to a bin size resolution  

  Int_t bins[dim] = {res, res, res};

  std::cout << "Bins defined \n" << std::endl;

  Double_t maxs[dim] = { 3500.,  3500., 3558.2+3500.}; // Detector size definitions
  Double_t mins[dim] = {-3500., -3500., 3558.2-3500.}; // 7000 represents 7m

  THnSparseF* voxels = new THnSparseF("Voxels", "", dim, bins, mins, maxs);
  VoxelsTree->Branch("voxels", "THnSparseF", &voxels);

  std::cout << "Voxels created and put into new branch \n" << std::endl; */

/*//////////////////
  for(Int_t i=0; i<TrackingTree->GetEntries(); ++i){

    TrackingTree->GetEntry(i); 

    voxels->Reset();

    // TPC Hits

    HitCollection tpcSdHits = SimulData->getTpcFidHits(); 
    
    for (Int_t j=0; j<tpcSdHits.size(); ++j) {
      SDHit tmpHit = tpcSdHits.at(j);
      Double_t xyz[3] = {tmpHit.getPosition().X(), tmpHit.getPosition().Y(), tmpHit.getPosition().Z()};
      voxels->Fill(xyz, tmpHit.getEdep()*1000.);
    }
    
    VoxelsTree->Fill();
  }
*//////////////////  

 
  std::cout << "For loops completed" << std::endl;
  
  VoxelsTree->Print();
  
  std::cout << "VoxelTree printed succesfully\n" << std::endl;  

//  TFile f1("~/Documents/TRex/hptrex/voxelresults.root", "NEW");
 f1.Write();
//  f1 = VoxelsTree->GetCurrentFile();
//  VoxelsTree->Write();

  f1.Close();  //PAULA: Attempt made to define output file and write to it
  return 0;
}
