#include "TCCDSimLoader.hxx"
#include "TKey.h"

trex::TCCDSimLoader::TCCDSimLoader(std::string inputFile){
  
  fFile=new TFile(inputFile.c_str());
  
  Detector = new TH3D("Detector", "Detector", 350, -350, 350, 350, -350, 350, 350, -350, 350);
  Detector->SetDirectory(0);

  TIter next(fFile->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())) {
    if(strstr(key->GetName(),"hReadOutIdeal")){
      TH2D* evt=(TH2D*)fFile->Get(key->GetName());
      fEvents.push_back(evt);
    }
  }
}


void trex::TCCDSimLoader::LoadEvent(unsigned int i){
  
  std::cout << "LOADING EVENT # " << i << std::endl;

  //Delete hits from last event
  for(std::vector<trex::TTPCHitPad*>::iterator hitPadIter=fHits.begin();hitPadIter!=fHits.end();++hitPadIter){
    delete *hitPadIter;
  }
  
  //Delete voxels from last event
  for(std::vector<voxel*>::iterator voxelIter=fVoxels.begin(); voxelIter!=fVoxels.end();++voxelIter){
    delete *voxelIter;
  }
      
  fHits.clear();
  fVoxels.clear();
  
  std::cout << " PEVIOUS EVENT HAS BEEN CLEANED UP " << std::endl;

  TH2D* evt=fEvents[i];
  
  for(int iy=0;iy<evt->GetNbinsX();++iy){
    for(int iz=0;iz<evt->GetNbinsY();++iz){
      float binContent=(evt->GetBinContent(iy,iz)+evt->GetBinContent(iy,iz+1)+evt->GetBinContent(iy,iz-1))/3.;
      if(binContent>1e-7){
	trex::TCCDSimLoader::voxel* voxel=new trex::TCCDSimLoader::voxel;
	voxel->x_pos=evt->GetXaxis()->GetBinCenter(iy)/2.;
	voxel->y_pos=evt->GetYaxis()->GetBinCenter(iz)/2.;
	voxel->z_pos=0.;
	voxel->Edep=binContent;
	voxel->time=0.;

	fVoxels.push_back(voxel);

	TLorentzVector pos4((*voxel).x_pos, (*voxel).y_pos, (*voxel).z_pos, (*voxel).time);
	trex::TTPCHitPad* hitPadPtr=new TTPCHitPad((*voxel).Edep,pos4);
	fHits.push_back(hitPadPtr);

	Detector->Fill(pos4.X(), pos4.Y(), pos4.Z());
      }
    }
  }  
}
