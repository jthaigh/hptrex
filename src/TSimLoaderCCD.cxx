#include "TSimLoaderCCD.hxx"


trex::TSimLoader::TSimLoader(std::string inputFile){

  Detector = new TH3D("Detector", "Detector", 1529, -710., 710., 1529, -710., 710., 1, 0, 1);
  Detector->SetDirectory(0);
  
  fFile=new TFile(inputFile.c_str(), "UPDATE");
  
  fVoxelsTree=(TTree*)fFile->Get("VoxelsTree");
  fBP_voxels=0;
  fBP_truthmatch_voxels=0;
  fVoxelsTree->SetBranchAddress("voxels", &fBP_voxels);
  fVoxelsTree->SetBranchAddress("truthmatch_voxels", &fBP_truthmatch_voxels);

  fTruthTree=(TTree*)fFile->Get("TruthTree");
  fTruthTree->SetBranchAddress("Momentum",fBP_Momentum);
  fTruthTree->SetBranchAddress("Xi",fBP_Xi);
  fTruthTree->SetBranchAddress("Yi",fBP_Yi);
  fTruthTree->SetBranchAddress("Zi",fBP_Zi);
  fTruthTree->SetBranchAddress("Xf",fBP_Xf);
  fTruthTree->SetBranchAddress("Yf",fBP_Yf);
  fTruthTree->SetBranchAddress("Zf",fBP_Zf);
  fTruthTree->SetBranchAddress("pdg",fBP_pdg);
  fTruthTree->SetBranchAddress("TrackID",fBP_TrackID);
  fTruthTree->SetBranchAddress("ParentID",fBP_ParentID);
  fTruthTree->SetBranchAddress("NParticles",&fBP_NParticles);
  
}


void trex::TSimLoader::LoadEvent(unsigned int i){

  /////////////////////////////////////
  //Delete everything from last event//
  /////////////////////////////////////

  for(std::vector<trex::TTPCHitPad*>::iterator hitPadIter=fHits.begin();hitPadIter!=fHits.end();++hitPadIter){
    delete *hitPadIter;
  }
  
  //Delete voxels from last event
  for(std::vector<voxel*>::iterator voxelIter=fVoxels.begin(); voxelIter!=fVoxels.end();++voxelIter){
    delete *voxelIter;
  }
  
  
  //Delete TrueHits from last event
  for(std::vector<TTrueHit*>::iterator TrueHitIter=fTrueHits.begin(); TrueHitIter!=fTrueHits.end();++TrueHitIter){
    delete *TrueHitIter;
  }
  
  //Delete TrueTracks from last event
  for(std::vector<trex::TTrueTrack*>::iterator truetrackIter=fTrueTracks.begin(); truetrackIter!=fTrueTracks.end(); ++truetrackIter){
    delete *truetrackIter;
  }

  fHits.clear();
  fVoxels.clear();
  fTrueHits.clear();
  fTrueTracks.clear();



  //////////////////////////////////////////////
  //Load event from trees into branch pointers//
  //////////////////////////////////////////////

  std::cout << "LOADING SPILL # " << i << std::endl;

  fVoxelsTree->GetEntry(i);
  fTruthTree->GetEntry(i);  



  //////////////////////////
  //True track information//
  //////////////////////////

  for(unsigned int iPart=0;iPart<fBP_NParticles;++iPart){
    fTrueTracks.push_back(new TTrueTrack());
    
    fTrueTracks.back()->SetEntries(fBP_pdg[iPart],
				   iPart,
				   fBP_TrackID[iPart],
				   fBP_ParentID[iPart], 
				   TVector3(fBP_Xi[iPart],fBP_Yi[iPart],fBP_Zi[iPart]),
				   TVector3(fBP_Xi[iPart],fBP_Yi[iPart],fBP_Zi[iPart]),
				   fBP_Momentum[iPart], 
				   fBP_NParticles);
  }

  std::vector<int> nHitsByTrack(fBP_NParticles,0);

  /////////////////////////////////////
  //Build voxels from input histogram//
  /////////////////////////////////////

  Int_t nVoxels = fBP_voxels->GetNbins(); 
  
  std::cout << "There are " << nVoxels << " Voxels in this Event" << std::endl;
  std::cout << "Of which " << fBP_truthmatch_voxels->GetNbins() << " are truth-matched" << std::endl;

  //JTH: THnSparse::GetBinContent(Int_t*) appears to be broken so I am loading the values into a map.
  std::map<int,short> truthMap;

  for(int i=0;i<fBP_truthmatch_voxels->GetNbins();++i){
    Int_t coords[3];
    Short_t val = fBP_truthmatch_voxels->GetBinContent(i,coords);
    truthMap[coords[0]*2000+coords[1]]=val;
  }

  for(int linInd=0; linInd<nVoxels; ++linInd){
          
    fVoxels.push_back(new voxel);
    voxel * voxelPtr = fVoxels.back();

    //Extracting information from the THnSparseF
   
    //Hard-coded Histogram dimensions (this is not good)
    Int_t coords[3];
    Double_t position[3];
        
    //PRD original resolution
    //double res = 2.34;
    //Int_t bins[3] = {513, 513, 1};
    //Double_t maxs[3] = { 600.21,  600.21, 0.};
    //Double_t mins[3] = {-600.21, -600.21, 1.};

    (*voxelPtr).Edep = fBP_voxels->GetBinContent(linInd, coords);
    (*voxelPtr).time = 0; //setting time to 0 for now until we have a T0 from other subdetectors

    Int_t voxelTrueTrack = (truthMap.count(coords[0]*2000+coords[1]))?(truthMap[coords[0]*2000+coords[1]]-1):-1;
    
    //Translate coordinates into positions
    //This requires a bit more thought. For now position = coordinates works. 
    for(int dim=0; dim<3; ++dim){
      //position[dim]=(mins[dim]) + res*coords[dim]; //real position (not good for cell-ID map right now)
      position[dim]=coords[dim];
    }

    //Swapping coordinates here because t2k beam direction is z and drift direction is x, whereas PRD beam direction is y and drift direction is z. 
    (*voxelPtr).x_pos = position[2];
    (*voxelPtr).y_pos = position[0];
    (*voxelPtr).z_pos = position[1];

    //std::cout << "Have found Voxel at position: " << coords[0] << " : " << coords[1] << " : " << coords[2] << std::endl;
    
    TLorentzVector pos4((*voxelPtr).x_pos, (*voxelPtr).y_pos, (*voxelPtr).z_pos, (*voxelPtr).time);

    trex::TTPCHitPad* hitPadPtr=new TTPCHitPad((*voxelPtr).Edep,pos4);
    if(voxelTrueTrack>=0){
      hitPadPtr->SetTrueTrack(fTrueTracks[voxelTrueTrack]);
      (nHitsByTrack[voxelTrueTrack])++;
    }

    fHits.push_back(hitPadPtr);

    Detector->Fill(pos4.X(), pos4.Y(), pos4.Z());

  }

  for(unsigned int iTrack=0;iTrack<fBP_NParticles;++iTrack){
    fTrueTracks[iTrack]->SetNumberOfHits(nHitsByTrack[iTrack]);
    std::cout<<"Track "<<iTrack<<" has "<<nHitsByTrack[iTrack]<<" hits."<<std::endl;
  }
 
}


unsigned int trex::TSimLoader::GetNEvents(){
  return fVoxelsTree->GetEntries();
}
