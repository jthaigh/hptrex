// eddy
#include "TTPCPathVolume.hxx"

ND::TTPCPathVolume::TTPCPathVolume(ND::TTPCUnitVolume* unitVolume){
  fUnitVolume = unitVolume;
  fFriends = std::vector<ND::TTPCUnitVolume*>();

  fXMin = 0;
  fXMax = 0;
  fXSize = 0;
  fYMin = 0;
  fYMax = 0;
  fYSize = 0;
  fZMin = 0;
  fZMax = 0;
  fZSize = 0;
  fAverageTime = 0.;
  fAveragePos = TVector3();
  fAveragePosXYZ = TVector3();

  fClosed = false;
  fIsXCluster = false;
}
ND::TTPCPathVolume::~TTPCPathVolume(){
}

std::vector<ND::TTPCUnitVolume*> ND::TTPCPathVolume::GetExtendedCell(int filter){
  // list to return
  std::vector<ND::TTPCUnitVolume*> extendedCells = std::vector<ND::TTPCUnitVolume*>();
  // add friends to the list
  if(GetHasCluster()){
    if((fIsVertical && (filter==0 || filter==1)) || (!fIsVertical && (filter==0 || filter==2))){
      for(std::vector<ND::TTPCUnitVolume*>::iterator ffriend = fFriends.begin(); ffriend != fFriends.end(); ++ffriend) extendedCells.push_back(*ffriend);
    }
  };

  return extendedCells;
}
ND::THandle<ND::TTPCHVCluster> ND::TTPCPathVolume::GetHits(){
  // cells whos hits to add
  std::vector<ND::TTPCUnitVolume*> volsToAdd = GetExtendedCell();
  ND::THitSelection hitsel;

  for(std::vector<ND::TTPCUnitVolume*>::iterator volIt = volsToAdd.begin(); volIt != volsToAdd.end(); ++volIt){
    ND::TTPCUnitVolume* vol = *volIt;
    if(!vol) continue;

    std::vector< ND::THandle<ND::TTPCHitPad> > hits;
    hits = vol->GetHits();
    for(std::vector< ND::THandle<ND::TTPCHitPad> >::iterator hitIt = vol->GetHitsBegin(); hitIt != vol->GetHitsEnd(); ++hitIt){
      ND::THandle<ND::TTPCHitPad> nhit = *hitIt;
      if (nhit)
        hitsel.AddHit(nhit);
      else
        throw;
    };
  };

  // null
  if(!hitsel.size()) return ND::THandle<ND::TTPCHVCluster>();
  // otherwise
  return ND::THandle<ND::TTPCHVCluster> ( new ND::TTPCHVCluster(hitsel, fIsVertical) );
}

void ND::TTPCPathVolume::MarkFriend(std::vector<ND::TTPCUnitVolume*>::iterator focusFriendIt){
  fClosed = false;
  *focusFriendIt = 0;
}
void ND::TTPCPathVolume::ClearMarked(){
  std::vector<ND::TTPCUnitVolume*>::iterator reaper;
  reaper = fFriends.begin();
  while(reaper != fFriends.end()){
    if(!(*reaper)) fFriends.erase(reaper);
    else(reaper++);
  };
}
void ND::TTPCPathVolume::ClearFriends(){
  fFriends.empty();
}

TVector3 ND::TTPCPathVolume::GetAvgPosXYZ(){
  TVector3 avgPosXYZ (0., 0., 0.);
  int norm = 0;

  for(std::vector<ND::TTPCUnitVolume*>::iterator friendVolIt = fFriends.begin(); friendVolIt != fFriends.end(); ++friendVolIt){
    ND::TTPCUnitVolume* friendVol = *friendVolIt;

    int weight = friendVol->size();
    norm += weight;
    avgPosXYZ += weight * friendVol->GetPosXYZ();
  };
  avgPosXYZ *= 1./(float)norm;

  return avgPosXYZ;
}

void ND::TTPCPathVolume::Close(){
  if(fClosed) return;

  fXMin = +999999;
  fXMax = -999999;
  fYMin = +999999;
  fYMax = -999999;
  fZMin = +999999;
  fZMax = -999999;

  fAverageTime = 0.;
  fAveragePos = TVector3(0., 0., 0.);
  fAveragePosXYZ = TVector3(0., 0., 0.);
  int norm = 0;
  for(std::vector<ND::TTPCUnitVolume*>::iterator friendIt = fFriends.begin(); friendIt != fFriends.end(); ++friendIt){
    ND::TTPCUnitVolume* frnd = *friendIt;

    fXMin = std::min(fXMin, frnd->GetX());
    fXMax = std::max(fXMax, frnd->GetX());
    fYMin = std::min(fYMin, frnd->GetY());
    fYMax = std::max(fYMax, frnd->GetY());
    fZMin = std::min(fZMin, frnd->GetZ());
    fZMax = std::max(fZMax, frnd->GetZ());

    fAverageTime += frnd->GetTime();
    fAveragePos += frnd->GetPos();
    fAveragePosXYZ += frnd->GetPosXYZ();
    norm ++;
  };
  fXSize = (fXMax - fXMin) + 1;
  fYSize = (fYMax - fYMin) + 1;
  fZSize = (fZMax - fZMin) + 1;

  fAverageTime /= double(norm);
  fAveragePos *= 1. / double(norm);
  fAveragePosXYZ *= 1. / double(norm);

  fClosed = true;
}

void ND::TTPCPathVolume::PrintPositions(bool showOrientation){
  std::cout << "( ";
  for(std::vector<ND::TTPCUnitVolume*>::iterator friendIt = fFriends.begin(); friendIt != fFriends.end(); ++friendIt){
    ND::TTPCUnitVolume* frnd = *friendIt;
    std::cout << frnd->GetX() << "," << frnd->GetY() << "," << frnd->GetZ() << " ";
  };
  std::cout << ")";
  if(showOrientation){
    std::cout << (fIsVertical ? "v" : "h");
    if(fIsXCluster) std::cout << "x";
  };
  std::cout << std::endl;
}
void ND::TTPCPathVolume::CheckHits(){
  // cells whos hits to check
  std::vector<ND::TTPCUnitVolume*> volsToCheck = GetExtendedCell();

  for(std::vector<ND::TTPCUnitVolume*>::iterator volIt = volsToCheck.begin(); volIt != volsToCheck.end(); ++volIt){
    ND::TTPCUnitVolume* vol = *volIt;
    if(!vol){
      if(ND::tpcDebug().PatternRecognition(DB_ERROR)) std::cout << "WARNING:  empty ND::TTPCUnitVolume* found in TTPCPathVolume" << std::endl;
    }
    else{
      if(ND::tpcDebug().PatternRecognition(DB_VERBOSE)) std::cout << "Attempting to access ND::TTPCUnitVolume hits…" << std::endl;
      std::vector< ND::THandle<ND::TTPCHitPad> > hits = vol->GetHits();
      if(ND::tpcDebug().PatternRecognition(DB_VERBOSE)) std::cout << "…success!" << std::endl;
      for(std::vector< ND::THandle<ND::TTPCHitPad> >::iterator hitIt = vol->GetHitsBegin(); hitIt != vol->GetHitsEnd(); ++hitIt){
        ND::THandle<ND::TTPCHitPad> nhit = *hitIt;
        if (!nhit){
          if(ND::tpcDebug().PatternRecognition(DB_ERROR)) std::cout << "WARNING:  non ND::TTPCHitPad hit found in TTPCPathVolume" << std::endl;
        };
      };
    };
  };
}
