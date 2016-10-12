// eddy
#include "TTPCPathVolume.hxx"

trex::TTPCPathVolume::TTPCPathVolume(trex::TTPCUnitVolume* unitVolume){
  fUnitVolume = unitVolume;

  fXMin = 0;
  fXMax = 0;
  fXSize = 0;
  fYMin = 0;
  fYMax = 0;
  fYSize = 0;
  fZMin = 0;
  fZMax = 0;
  fZSize = 0;

  fClosed = false;
  fIsXCluster = false;
}
trex::TTPCPathVolume::~TTPCPathVolume(){

}

std::vector<trex::TTPCUnitVolume*> trex::TTPCPathVolume::GetExtendedCell(int filter){
  // list to return
  std::vector<trex::TTPCUnitVolume*> extendedCells = std::vector<trex::TTPCUnitVolume*>();
  // add friends to the list
  if(GetHasCluster()){
    if((fIsVertical && (filter==0 || filter==1)) || (!fIsVertical && (filter==0 || filter==2))){
      for(std::vector<trex::TTPCUnitVolume*>::iterator ffriend = fFriends.begin(); ffriend != fFriends.end(); ++ffriend) extendedCells.push_back(*ffriend);
    }
  };

  return extendedCells;
}

std::vector<trex::TTPCHitPad*> trex::TTPCPathVolume::GetHits(){
  
  std::vector<trex::TTPCHitPad*> output;
  // cells whos hits to add
  std::vector<trex::TTPCUnitVolume*> volsToAdd = GetExtendedCell();
  
  for(std::vector<trex::TTPCUnitVolume*>::iterator volIt = volsToAdd.begin(); volIt != volsToAdd.end(); ++volIt){
    trex::TTPCUnitVolume* vol = *volIt;
    if(!vol) continue;

    std::vector< trex::TTPCHitPad* > hits;
    hits = vol->GetHits();
    for(std::vector< trex::TTPCHitPad* >::iterator hitIt = vol->GetHitsBegin(); hitIt != vol->GetHitsEnd(); ++hitIt){
      trex::TTPCHitPad* nhit = *hitIt;
      if (nhit)
        output.push_back(nhit);
      else
        throw;
    }
  }

  return output;
}

void trex::TTPCPathVolume::MarkFriend(std::vector<trex::TTPCUnitVolume*>::iterator focusFriendIt){
  fClosed = false;
  *focusFriendIt = 0;
}
void trex::TTPCPathVolume::ClearMarked(){
  std::vector<trex::TTPCUnitVolume*>::iterator reaper;
  reaper = fFriends.begin();
  while(reaper != fFriends.end()){
    if(!(*reaper)) fFriends.erase(reaper);
    else(reaper++);
  };
}
void trex::TTPCPathVolume::ClearFriends(){
  fFriends.clear();
}

TVector3 trex::TTPCPathVolume::GetAvgPosXYZ(){
  TVector3 avgPosXYZ (0., 0., 0.);
  int norm = 0;

  for(std::vector<trex::TTPCUnitVolume*>::iterator friendVolIt = fFriends.begin(); friendVolIt != fFriends.end(); ++friendVolIt){
    trex::TTPCUnitVolume* friendVol = *friendVolIt;

    int weight = friendVol->size();
    norm += weight;
    avgPosXYZ += weight * friendVol->GetPosXYZ();
  };
  avgPosXYZ *= 1./(float)norm;

  return avgPosXYZ;
}

void trex::TTPCPathVolume::Close(){
  if(fClosed) return;

  fXMin = +999999;
  fXMax = -999999;
  fYMin = +999999;
  fYMax = -999999;
  fZMin = +999999;
  fZMax = -999999;

  fAveragePos = TVector3(0., 0., 0.);
  fAveragePosXYZ = TVector3(0., 0., 0.);
  int norm = 0;
  for(std::vector<trex::TTPCUnitVolume*>::iterator friendIt = fFriends.begin(); friendIt != fFriends.end(); ++friendIt){
    trex::TTPCUnitVolume* frnd = *friendIt;

    fXMin = std::min(fXMin, frnd->GetX());
    fXMax = std::max(fXMax, frnd->GetX());
    fYMin = std::min(fYMin, frnd->GetY());
    fYMax = std::max(fYMax, frnd->GetY());
    fZMin = std::min(fZMin, frnd->GetZ());
    fZMax = std::max(fZMax, frnd->GetZ());

    fAveragePos += frnd->GetPos();
    fAveragePosXYZ += frnd->GetPosXYZ();
    norm ++;
  };
  fXSize = (fXMax - fXMin) + 1;
  fYSize = (fYMax - fYMin) + 1;
  fZSize = (fZMax - fZMin) + 1;

  fAveragePos *= 1. / double(norm);
  fAveragePosXYZ *= 1. / double(norm);

  fClosed = true;
}

void trex::TTPCPathVolume::PrintPositions(bool showOrientation){
  std::cout << "( ";
  for(std::vector<trex::TTPCUnitVolume*>::iterator friendIt = fFriends.begin(); friendIt != fFriends.end(); ++friendIt){
    trex::TTPCUnitVolume* frnd = *friendIt;
    std::cout << frnd->GetX() << "," << frnd->GetY() << "," << frnd->GetZ() << " ";
  };
  std::cout << ")";
  if(showOrientation){
    std::cout << (fIsVertical ? "v" : "h");
    if(fIsXCluster) std::cout << "x";
  };
  std::cout << std::endl;
}
void trex::TTPCPathVolume::CheckHits(){
  // cells whos hits to check
  std::vector<trex::TTPCUnitVolume*> volsToCheck = GetExtendedCell();

  for(std::vector<trex::TTPCUnitVolume*>::iterator volIt = volsToCheck.begin(); volIt != volsToCheck.end(); ++volIt){
    trex::TTPCUnitVolume* vol = *volIt;
    if(!vol){
      std::cout << "WARNING:  empty trex::TTPCUnitVolume* found in TTPCPathVolume" << std::endl;
    }
    else{
      std::cout << "Attempting to access trex::TTPCUnitVolume hits…" << std::endl;
      std::vector< trex::TTPCHitPad* > hits = vol->GetHits();
      std::cout << "…success!" << std::endl;
      for(std::vector< trex::TTPCHitPad* >::iterator hitIt = vol->GetHitsBegin(); hitIt != vol->GetHitsEnd(); ++hitIt){
        trex::TTPCHitPad* nhit = *hitIt;
        if (!nhit){
          std::cout << "WARNING:  non trex::TTPCHitPad hit found in TTPCPathVolume" << std::endl;
        }
      }
    }
  }
}
