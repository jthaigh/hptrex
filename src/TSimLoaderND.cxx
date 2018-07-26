#include "TSimLoaderND.hxx"


trex::TSimLoaderND::TSimLoaderND(std::string inputFile){
  
  fFile=new TFile(inputFile.c_str(), "UPDATE");
  
  fTree=(TTree*)fFile->Get("EventRecord");
  //fTree->SetBranchStatus("event",0); //this just needs to be there until the genie libs are linked

  //fVoxelsTree=(TTree*)fFile->Get("VoxelsTree");

  //fReconTree=new TTree("TPCRecon", "TPCRecon");

   
  fEventRecordBranch=0;
  fTree->SetBranchAddress("EventRecord",&fEventRecordBranch);

  //fTree->SetBranchAddress("GeantTrackingTruth", &fGeantBranch);

  //fVoxelBranch=0;
  //fVoxelsTree->SetBranchAddress("voxels", &fVoxelBranch);

  Detector = new TH3D("Detector", "Detector", 350, -350, 350, 350, -350, 350, 350, 355.82-350, 355.82+350);
  Detector->SetDirectory(0);
}


void trex::TSimLoaderND::LoadEvent(unsigned int i){
  
  std::cout << "LOADING EVENT # " << i << std::endl;
  
  //Delete hits from last even
  for(std::vector<trex::TTPCHitPad*>::iterator hitPadIter=fHits.begin();hitPadIter!=fHits.end();++hitPadIter){
    delete *hitPadIter;
  }
  
  //Delete voxels from last event
  for(std::vector<voxel*>::iterator voxelIter=fVoxels.begin(); voxelIter!=fVoxels.end();++voxelIter){
    delete *voxelIter;
  }

  //Delete TrueTracks from last event
  for(auto TrueTrackIter=fTrueTracks.begin(); TrueTrackIter!=fTrueTracks.end();++TrueTrackIter){
    delete *TrueTrackIter;
  }  
  
  //Delete TrueHits from last event
  for(std::vector<TTrueHit*>::iterator TrueHitIter=fTrueHits.begin(); TrueHitIter!=fTrueHits.end();++TrueHitIter){
    delete *TrueHitIter;
  }
  
  fHits.clear();
  fVoxels.clear();
  fTrueTracks.clear();
  fTrueHits.clear();

  //fVoxelsTree->GetEntry(i);
  fTree->GetEntry(i);

  const std::vector<gastpc::MCParticle*>& MCParticles = fEventRecordBranch->GetMCParticles();
  const std::vector<gastpc::MCTrack*>& MCTracks = fEventRecordBranch->GetMCTracks();

  //HitCollection simHits = fSimulDataBranch->getTpcFidHits();
    
  Int_t nEntries = MCTracks.size(); 
  
  std::cout << "This Event contains " << nEntries << " Tracks " << std::endl;
  
  for(auto trackIt=MCTracks.begin(); trackIt!=MCTracks.end(); ++trackIt){
    
    const gastpc::MCTrack* theTrack=*trackIt;

    if(theTrack->GetLabel()!=std::string("TPC")) continue;

    gastpc::MCParticle* theParticle=theTrack->GetMCParticle();
    gastpc::Vector4D initPos=theParticle->GetInitialXYZT();
    gastpc::Vector4D finalPos=theParticle->GetFinalXYZT();
    gastpc::Vector3D initMom=theParticle->GetInitialMomentum(); 
    int mcID=theParticle->GetMCID();
    double initMom_mag=TMath::Sqrt(initMom.GetX()*initMom.GetX()+
			    initMom.GetY()*initMom.GetY()+
			    initMom.GetZ()*initMom.GetZ());
    

    fTrueTracks.emplace_back();
    fTrueTracks.back()=new trex::TTrueTrack;
    fTrueTracks.back()->SetEntries(theParticle->GetPDGCode(),
				   mcID,
				   mcID,
				   /*0,*/
				   (theParticle->GetParent()?theParticle->GetParent()->GetMCID():-1),
				   TVector3(initPos.GetX(),initPos.GetY(),initPos.GetZ()),
				   TVector3(finalPos.GetX(),finalPos.GetY(),finalPos.GetZ()),
				   initMom_mag,0);


    const std::vector<gastpc::MCHit*>& MCHits = (*trackIt)->GetMCHits();

    for(auto hitIt=MCHits.begin(); hitIt!=MCHits.end(); ++hitIt){
      
      gastpc::Vector4D xyzt = (*hitIt)->GetXYZT();
      double amplitude = (*hitIt)->GetAmplitude();

      if(xyzt.GetX()<-1240||xyzt.GetX()>1240||
	 xyzt.GetY()<-1240||xyzt.GetY()>1240||
	 xyzt.GetZ()<550||xyzt.GetZ()>7060){
	continue;
      }

      //SDHit& hit=*hitIter; does not make sense here right now
      //double hitTime=hit.getPosition().T();
      
      fVoxels.push_back(new voxel);
      voxel * voxelPtr = fVoxels.back();
      
      //Extracting information from the THnSparseF
      
    
      //Hard-coded Histogram dimensions (this is not good)
      //Int_t coords[3];
      //Double_t position[3];
      
      //DUNE ND
      //Int_t bins[3] = {700, 700, 700};
      //Double_t maxs[3] = { 3500.,  3500., 3558.2+3500.};
      //Double_t mins[3] = {-3500., -3500., 3558.2-3500.};
            
      (*voxelPtr).Edep = amplitude;
      (*voxelPtr).time = 0; //setting time to 0 for now until we have a T0 from other subdetectors
      
      /*
      //Translate coordinates into positions
      for(int dim=0; dim<3; ++dim){
	//position[dim]=(mins[dim]) + res*coords[dim]; //real position (not good for cell-ID map right now
	position[dim]=coords[dim];
      }
      */
      
      (*voxelPtr).x_pos = xyzt.GetX()*0.1;
      (*voxelPtr).y_pos = xyzt.GetY()*0.1;
      (*voxelPtr).z_pos = xyzt.GetZ()*0.1;
      
      //std::cout << "Have found Voxel at position: " << (*voxelPtr).x_pos << " : " << (*voxelPtr).y_pos << " : " << (*voxelPtr).z_pos << std::endl;
      
      TLorentzVector pos4((*voxelPtr).x_pos, (*voxelPtr).y_pos, (*voxelPtr).z_pos, (*voxelPtr).time);
      
      trex::TTPCHitPad* hitPadPtr=new TTPCHitPad((*voxelPtr).Edep,pos4);
      hitPadPtr->SetTrueTrack(fTrueTracks.back());
      
      fHits.push_back(hitPadPtr);
      
      Detector->Fill(pos4.X(), pos4.Y(), pos4.Z());
      
    }
  }
  
  //fTree->GetEntry(i);

  /*HitCollection simHits = fSimulDataBranch->getTpcFidHits();
  
  for(HitCollection::iterator hitIter=simHits.begin();hitIter!=simHits.end();++hitIter){

    //std::cout << "Entering the HitCollection LOOP!" << std::endl;x
    
  SDHit& hit=*hitIter;
  double TrueEdep=hit.getEdep();
  TLorentzVector TruePos4=hit.getPosition();
  int TrueTrackID = hit.getTrackID();
  //std::cout << "TRACK ID: " << TrueTrackID << std::endl;
  int pdg = hit.getPDG();
  int charge = hit.getCharge();
    
  fTrueHits.push_back(new TTrueHit());
  TTrueHit * TrueHitPtr = fTrueHits.back();
  (*TrueHitPtr).TrueEdep = TrueEdep;
  (*TrueHitPtr).TruePos4 = TruePos4;
  (*TrueHitPtr).pdg = pdg;
  (*TrueHitPtr).TrueTrackID = TrueTrackID;
  (*TrueHitPtr).charge = charge;
  //delete TrueHitPtr;
  }
  
  std::sort(fTrueHits.begin(), fTrueHits.end());
  */
  
}


unsigned int trex::TSimLoaderND::GetNEvents(){
  return fTree->GetEntries();
}




