#include "TSimLoaderCCD.hxx"


trex::TSimLoader::TSimLoader(std::string inputFile){
  
  fFile=new TFile(inputFile.c_str(), "UPDATE");
  

  fVoxelsTree=(TTree*)fFile->Get("VoxelsTree");

  fVoxelBranch=0;
  fVoxelsTree->SetBranchAddress("voxels", &fVoxelBranch);

  Detector = new TH3D("Detector", "Detector", 300, -600.21, 600.21, 350, -600.21, 600.21, 1, 0, 1);
  Detector->SetDirectory(0);
}


void trex::TSimLoader::LoadEvent(unsigned int i){

  Int_t ImageNumber=0;
  fVoxelsTree->SetBranchAddress("ImageNumber", &ImageNumber);

  Double_t Momentum=0;
  fVoxelsTree->SetBranchAddress("Momentum", &Momentum);

  TVector3* TrueXi=0;
  fVoxelsTree->SetBranchAddress("Xi", &TrueXi);

  TVector3* TrueXf=0;
  fVoxelsTree->SetBranchAddress("Xf", &TrueXf);
  
  Int_t PDG=0;
  fVoxelsTree->SetBranchAddress("PDG", &PDG);
  
  Int_t TrackID = 0;
  fVoxelsTree->SetBranchAddress("TrackID", &TrackID);
  
  Int_t ParentID = 0;
  fVoxelsTree->SetBranchAddress("ParentID", &ParentID);
  
  Int_t ProOrPi = 0;
  fVoxelsTree->SetBranchAddress("ProOrPi", &ProOrPi);


  std::cout << "LOADING EVENT # " << i << std::endl;

  fVoxelsTree->GetEntry(i);

  std::cout << "Processing Image # " << ImageNumber << std::endl;
  
   
  //Clean up previous event but only if we are actually processing a new image - if not just add another single track to the event.
  if(ImageNumber != fimageNumber){

    //Delete hits from last event
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

    fimageNumber = ImageNumber;
    
  }

  fTrueTracks.push_back(new TTrueTrack());
 
  Int_t TrackNumber = fTrueTracks.size();
  std::cout << "This is Track number " << TrackNumber << " of Image number " << ImageNumber << std::endl;

  fTrueMultiplicity = TrackNumber;
  
  fTrueTracks.back()->SetEntries(PDG,TrackNumber,TrackID, ProOrPi, ParentID, *TrueXi, *TrueXf, Momentum);

  Int_t nVoxels = fVoxelBranch->GetNbins(); 
  
  std::cout << "There are " << nVoxels << " Voxels in this Event" << std::endl;

  if(TrackID==1){
    fTrueTracks.back()->SetNumberOfHits(nVoxels);}
  else{std::cout << "THIS IS A TRACK WITH HIGHER TRACK ID!!!" << std::endl;}
 
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

    (*voxelPtr).Edep = fVoxelBranch->GetBinContent(linInd, coords);
    (*voxelPtr).time = 0; //setting time to 0 for now until we have a T0 from other subdetectors

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
    hitPadPtr->SetTrueTrack(fTrueTracks.back());

    fHits.push_back(hitPadPtr);

    Detector->Fill(pos4.X(), pos4.Y(), pos4.Z());

  }

  //std::cout << "These Hits belong to a Track with TrackID: " << fTrueTracks.back()->GetTrackID() << std::endl;
 
}


unsigned int trex::TSimLoader::GetNEvents(){
  return fVoxelsTree->GetEntries();
}



