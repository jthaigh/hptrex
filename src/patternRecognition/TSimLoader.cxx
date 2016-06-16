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
    voxel * voxelPtr;

    //Extracting information from the THnSparseF
   
    Int_t coords[3];
    
    (*voxelPtr).Edep = fVoxelBranch->GetBinContent(linInd, coords);
    (*voxelPtr).time = 0; //setting time to 0 for now until we have a T0 from other subdetectors
    (*voxelPtr).x_pos = coords[0];
    (*voxelPtr).y_pos = coords[1];
    (*voxelPtr).z_pos = coords[2];

    std::cout << "Have found Voxel at position: " << coords[0] << " : " << coords[1] << " : " << coords[2] << std::endl;

    fVoxels.push_back(voxelPtr);
    
    TLorentzVector pos4((*voxelPtr).x_pos, (*voxelPtr).y_pos, (*voxelPtr).z_pos, (*voxelPtr).time);

    trex::TTPCHitPad* hitPadPtr=new TTPCHitPad((*voxelPtr).Edep,pos4);

    fHits.push_back(hitPadPtr);

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
