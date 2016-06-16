#include "TSimLoader.hxx"

trex::TSimLoader::TSimLoader(std::string inputFile){

  fFile=new TFile(inputFile.c_str());

  fTree=(TTree*)fFile->Get("Tracking");
  fVoxelsTree=(TTree*)fFile->Get("VoxelsTree");

  fSimulDataBranch=0;
  fTree->SetBranchAddress("SimulData",&fSimulDataBranch);

  fVoxelBranch=0;
  fVoxelsTree->SetBranchAddress("voxels", &fVoxelBranch);
}


void trex::TSimLoader::LoadEvent(unsigned int i){

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

  
  fVoxelsTree->GetEntry(i);

  std::cout << "JUST TRIED GETTING ENTRY FROM fTree" << std::endl;

  //HitCollection simHits = fSimulDataBranch->getTpcFidHits();
  //THnSparseF VoxelCollection = (THnSparseF) fVoxelBranch->Clone();   

  std::cout << "BRANCHES HAVE BEEN SET UP. ENTERING LOOP NOW." << std::endl;

  Int_t nVoxels = fVoxelBranch->GetNbins(); 

  std::cout << "FOUND " << nVoxels << " Voxels " << std::endl; 

  for(int linInd=0; linInd<nVoxels; ++linInd){

    //SDHit& hit=*hitIter;
    //double hitTime=hit.getPosition().T();
    
    fVoxels.push_back(new voxel);
    voxel * voxelPtr = fVoxels.back();

    //Extracting information from the THnSparseF
   
    Int_t coords[3];
    Double_t position[3];
    Int_t bins[3] = {700, 700, 700};
    Double_t maxs[3] = { 3500.,  3500., 3558.2+3500.};
    Double_t mins[3] = {-3500., -3500., 3558.2-3500.};
        
    (*voxelPtr).Edep = fVoxelBranch->GetBinContent(linInd, coords);
    (*voxelPtr).time = 0; //setting time to 0 for now until we have a T0 from other subdetectors

    //Translate coordinates into positions
    for(int dim=0; dim<3; ++dim){
      position[dim]=(mins[dim]/10) + coords[dim]; //position in cm
    }

    
    (*voxelPtr).x_pos = position[0];
    (*voxelPtr).y_pos = position[1];
    (*voxelPtr).z_pos = position[2];

    //std::cout << "Have found Voxel at position: " << coords[0] << " : " << coords[1] << " : " << coords[2] << std::endl;

    
    TLorentzVector pos4((*voxelPtr).x_pos, (*voxelPtr).y_pos, (*voxelPtr).z_pos, (*voxelPtr).time);

    trex::TTPCHitPad* hitPadPtr=new TTPCHitPad((*voxelPtr).Edep,pos4);

    fHits.push_back(hitPadPtr);

    fVoxels.at(linInd)->printVoxel(); 

  }

  //for(HitCollection::iterator hitIter=simHits.begin();hitIter!=simHits.end();++hitIter){
    //SDHit& hit=*hitIter;
    //double eDep=hit.getEdep();
    //TLorentzVector pos4=hit.getPosition();
    //trex::TTPCHitPad* hitPadPtr=new TTPCHitPad(eDep,pos4);
    //fHits.push_back(hitPadPtr);
    //}

}

unsigned int trex::TSimLoader::GetNEvents(){
  return fTree->GetEntries();
}
