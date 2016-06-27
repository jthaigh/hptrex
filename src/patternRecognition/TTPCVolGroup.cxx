// eddy
#include "TTPCVolGroup.hxx"

//ClassImp(trex::TTPCVolGroup);
trex::TTPCVolGroup::~TTPCVolGroup(){
}

unsigned int trex::TTPCVolGroup::sMaxID = 0;

trex::TTPCVolGroup::TTPCVolGroup(trex::TTPCLayout* layout, unsigned int id){

  fXLean = 0;
  fYLean = 0;
  fZLean = 0;

  fXMin = 9999 ;
  fXMax = -9999;
  fXSize = 0;
  fYMin = 9999 ;
  fYMax = -9999;
  fYSize = 0;
  fZMin = 9999 ;
  fZMax = -9999;
  fZSize = 0;

  fSigmaPadX = 0.;
  fSigmaPadY = 0.;
  fSigmaPadZ = 0.;
  fAveragePadID = 0;
  fCharge = 0.;
  fAverageCharge = 0.;
  fIsClosed = false;

  fLayout = layout;
  fID = id;
}

void trex::TTPCVolGroup::AddHitMap(std::map<long, trex::TTPCUnitVolume*> hitMap){
  // add map of hits to event
  fIsClosed = false;
  fHitMap = std::map<long, trex::TTPCUnitVolume*>(hitMap);
}
void trex::TTPCVolGroup::AddHits(std::map<long, trex::TTPCUnitVolume*> hitMap){
  // insert hits from input map to this group's map
  fIsClosed = false;
  fHitMap.insert(hitMap.begin(), hitMap.end());
}
void trex::TTPCVolGroup::AddHits(std::vector<trex::TTPCUnitVolume*> hitList){
  fIsClosed = false;
  AddHits(hitList.begin(), hitList.end());
}
void trex::TTPCVolGroup::AddHits(std::vector<trex::TTPCUnitVolume*>::iterator hitListBegin, std::vector<trex::TTPCUnitVolume*>::iterator hitListEnd){
  fIsClosed = false;
  for(std::vector<trex::TTPCUnitVolume*>::iterator volIt = hitListBegin; volIt != hitListEnd; ++volIt){
    trex::TTPCUnitVolume* vol = *volIt;
    fHitMap[vol->GetID()] = vol;
  };
}
void trex::TTPCVolGroup::AddHits(trex::TTPCVolGroup& hits){
  fIsClosed = false;
  // insert all hits in input group's map into this group's map
  for(std::map<long, trex::TTPCUnitVolume*>::iterator el = hits.begin(); el != hits.end(); ++el){
    fHitMap[el->first] = el->second;
  };
}
void trex::TTPCVolGroup::AddHits(std::map<long, trex::TTPCUnitVolume*>::iterator hitMapBegin, std::map<long, trex::TTPCUnitVolume*>::iterator hitMapEnd){
  fIsClosed = false;
  // insert all hits between iterators into this group's map
  for(std::map<long, trex::TTPCUnitVolume*>::iterator el = hitMapBegin; el != hitMapEnd; ++el){
    fHitMap[el->first] = el->second;
  };
}
bool trex::TTPCVolGroup::AddHit(trex::TTPCUnitVolume* hit, bool safe){
  fIsClosed = false;
  if (!hit) return false;
  // get hit id 
  long id = hit->GetID();
  // if it doesn't already exist in map, add it and return true
  if(!safe){
    fHitMap[id] = hit; 
    return true;
  }
  else if(!GetHit(id, true)){
    fHitMap[id] = hit;
    return true;
  };
  // otherwise, return false
  return false;
}

bool trex::TTPCVolGroup::RemoveHit(long id){
  fIsClosed = false;
  // try to find cell in current map
  std::map<long, trex::TTPCUnitVolume*>::iterator el = fHitMap.find(id);
  if (el == fHitMap.end()) return false;
  // if found, remove it
  fHitMap.erase(el);
  return true;
}

void trex::TTPCVolGroup::MarkHit(long id){
  fIsClosed = false;
  fHitMap[id] = 0;
}

void trex::TTPCVolGroup::ClearMarked(){
  fIsClosed = false;
  std::vector<long> keysToDelete = std::vector<long>();
  for(std::map<long, trex::TTPCUnitVolume*>::iterator hitEl = fHitMap.begin(); hitEl != fHitMap.end(); ++hitEl){
    if(!hitEl->second) keysToDelete.push_back(hitEl->first);
  };
  for(std::vector<long>::iterator keyIt = keysToDelete.begin(); keyIt != keysToDelete.end(); ++keyIt){
    fHitMap.erase(*keyIt);
  };
}

void trex::TTPCVolGroup::MergeHits(trex::TTPCVolGroup& hits){
  fIsClosed = false;
  for(std::map<long, trex::TTPCUnitVolume*>::iterator el = hits.begin(); el != hits.end(); ++el){
    std::map<long, trex::TTPCUnitVolume*>::iterator id = fHitMap.find(el->first);
    // if cell isn't found, add it to this map
    if (id == fHitMap.end()) fHitMap[el->first] = el->second;
    // if it is, add it's charge to the corresponding cell in this group's map
    else id->second->SetQ(id->second->GetQ() + el->second->GetQ());
  };
}

long trex::TTPCVolGroup::GetNearestHit(int x, int y, int z, int maxDist){
  int minDist = 999999;
  long id = -1;
  // loop over map looking for minimum distance
  for(std::map<long, trex::TTPCUnitVolume*>::iterator it = fHitMap.begin(); it != fHitMap.end(); ++it){
    int dX = it->second->GetX() - x;
    int dY = it->second->GetY() - y;
    int dZ = it->second->GetZ() - z;

    int dist = dX*dX + dY*dY + dZ*dZ; 
    if(dist > maxDist*maxDist) continue;
    if(dist > minDist) continue;
    minDist = dist;
    id = it->second->GetID();
  };
  // return correcponding cell id
  return id;
}

TVector3 trex::TTPCVolGroup::GetNearestHitPos(int x, int y, int z, int maxDist){
  // find nearest hit
  long id = GetNearestHit(x, y, z, maxDist);
  // if unfound return default vector
  if (id < 0) return TVector3(-9999., -9999., -9999.);
  // return its position
  return fHitMap[id]->GetPos(); 
}

std::vector<trex::TTPCHitPad*> trex::TTPCVolGroup::GetHits(){
  std::vector<trex::TTPCHitPad*> outHits;

  for(std::map<long, trex::TTPCUnitVolume*>::iterator el = fHitMap.begin(); el != fHitMap.end(); ++el){
    std::vector< trex::TTPCHitPad* > hits = el->second->GetHits();
    for(std::vector< trex::TTPCHitPad* >::iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt){
      trex::TTPCHitPad* hit = *hitIt;
      outHits.push_back(hit);
    }
  }
  return std::move(outHits);
}

bool trex::TTPCVolGroup::Contains(long id){
  // search this group for the provided id
  std::map<long, trex::TTPCUnitVolume*>::iterator el = fHitMap.find(id);
  return (el != fHitMap.end());
}

bool trex::TTPCVolGroup::Contains(TTPCUnitVolume* vol){
  return Contains(vol->GetID());
}

bool trex::TTPCVolGroup::GetLeanValid(trex::TTPCUnitVolume* vol){
  Close();
  if(fXLean != 0 && vol->GetX() != fAveragePad.x) return false;
  if(fYLean != 0 && vol->GetY() != fAveragePad.y) return false;
  if(fZLean != 0 && vol->GetZ() != fAveragePad.z) return false;
  return true;
}

trex::TTPCUnitVolume* trex::TTPCVolGroup::GetHit(long id, bool safe){
  safe = true;
  // try to find cell in current map
  std::map<long, trex::TTPCUnitVolume*>::iterator el;
  if(safe){
    // check if exists in hit map
    el = fHitMap.find(id);
    if (el == fHitMap.end()) return 0;
  };
  return fHitMap[id];
}

void trex::TTPCVolGroup::Close(){
  if(fIsClosed) return;
  if(fHitMap.size() < 1) return;

  fXMin = +9999;
  fXMax = -9999;
  fYMin = +9999;
  fYMax = -9999;
  fZMin = +9999;
  fZMax = -9999;
  for(std::map<long, trex::TTPCUnitVolume*>::iterator vol = fHitMap.begin(); vol != fHitMap.end(); ++vol){
    fXMin = std::min(fXMin, vol->second->GetX());
    fXMax = std::max(fXMax, vol->second->GetX());
    fYMin = std::min(fYMin, vol->second->GetY());
    fYMax = std::max(fYMax, vol->second->GetY());
    fZMin = std::min(fZMin, vol->second->GetZ());
    fZMax = std::max(fZMax, vol->second->GetZ());
  };

  fXSize = fXMax-fXMin + 1;
  fYSize = fYMax-fYMin + 1;
  fZSize = fZMax-fZMin + 1;

  // must be kept in this order (note to self, fix this at some point)
  FindAveragePadPosID();
  FindSigmaPads();
  FindCharge();

  fIsClosed = true;
}

void trex::TTPCVolGroup::FindSigmaPads(){
  trex::TTPCCellInfo3D avgPos = fAveragePad;
  // initialise sum of squared of differences to zero
  double sigmaSumX = 0.;
  double sigmaSumY = 0.;
  double sigmaSumZ = 0.;
  for(std::map<long, trex::TTPCUnitVolume*>::iterator it = fHitMap.begin(); it != fHitMap.end(); ++it){
    // get position difference in x, y or z
    int xDiff = avgPos.x - it->second->GetX();
    int yDiff = avgPos.y - it->second->GetY();
    int zDiff = avgPos.z - it->second->GetZ();

    // sum square of position difference
    sigmaSumX += xDiff*xDiff;
    sigmaSumY += yDiff*yDiff;
    sigmaSumZ += zDiff*zDiff;
  };
  // get square root and divide by size
  float sizeMult = 1./float(fHitMap.size());
  fSigmaPadX = std::sqrt( sizeMult * sigmaSumX );
  fSigmaPadY = std::sqrt( sizeMult * sigmaSumY );
  fSigmaPadZ = std::sqrt( sizeMult * sigmaSumZ );
}
void trex::TTPCVolGroup::FindAveragePadPosID(){
  // find averages in all co-ordinates
  int avgX = 0;
  int avgY = 0;
  int avgZ = 0;
  int avgN = 0;

  // iterate with progressively easier lean conditions (one iteration should usually be enough)
  for(int iter=4; iter > 0; iter--){
    for(std::map<long, trex::TTPCUnitVolume*>::iterator padEl = fHitMap.begin(); padEl != fHitMap.end(); ++padEl){
      // ignore if one of the lean filters is transgressed
      if(iter > 1){
        // try x first since most awkward conflicts come from spirals
        if(fXLean < 0 && padEl->second->GetX() > fXMin) continue;
        if(fXLean > 0 && padEl->second->GetX() < fXMax) continue;
      };
      if(iter > 2){
        if(fZLean < 0 && padEl->second->GetZ() > fZMin) continue;
        if(fZLean > 0 && padEl->second->GetZ() < fZMax) continue;
      };
      if(iter > 3){
        if(fYLean < 0 && padEl->second->GetY() > fYMin) continue;
        if(fYLean > 0 && padEl->second->GetY() < fYMax) continue;
      };

      avgX += padEl->second->GetX();
      avgY += padEl->second->GetY();
      avgZ += padEl->second->GetZ();
      avgN ++;
    };
    if(avgN > 0) break;
  };
  if(avgN == 0) return;

  avgX /= avgN;
  avgY /= avgN;
  avgZ /= avgN;

  // stick to pad in group
  double minDist2 = 99999999.;
  // iterate with progressively easier lean conditions (one iteration should usually be enough)
  for(int iter=4; iter > 0; iter--){
    for(std::map<long, trex::TTPCUnitVolume*>::iterator padEl = fHitMap.begin(); padEl != fHitMap.end(); ++padEl){
      if(iter > 1){
        if(fZLean < 0 && padEl->second->GetZ() > fZMin) continue;
        if(fZLean > 0 && padEl->second->GetZ() < fZMax) continue;
      };
      if(iter > 2){
        if(fYLean < 0 && padEl->second->GetY() > fYMin) continue;
        if(fYLean > 0 && padEl->second->GetY() < fYMax) continue;
      };
      if(iter > 3){
        if(fXLean < 0 && padEl->second->GetX() > fXMin) continue;
        if(fXLean > 0 && padEl->second->GetX() < fXMax) continue;
      };

      int dX = padEl->second->GetX() - avgX;
      int dY = padEl->second->GetY() - avgY;
      int dZ = padEl->second->GetZ() - avgZ;
      double curDist2 = dX*dX + dY*dY + dZ*dZ;

      if(curDist2 < minDist2){
        fAveragePadID = padEl->first;
        fAverageUnitVolume = padEl->second;
        fAveragePad = padEl->second->GetCellInfo3D();
        fAveragePosition = padEl->second->GetPos();
        minDist2 = curDist2;
      };
    };
    if(minDist2 < 9999.) break;
  };
}
void trex::TTPCVolGroup::FindCharge(){
  // initialise sum to zero
  float qSum = 0; 
  for(std::map<long, trex::TTPCUnitVolume*>::iterator el = fHitMap.begin(); el != fHitMap.end(); ++el){
    // add charge of every element in group
    qSum += el->second->GetQ();
  };
  fCharge = qSum;
  if(fHitMap.size() > 0){
    fAverageCharge = qSum / (float)fHitMap.size();
  };
}

//MDH
//Under new scheme (HitSel->vector of TTPCHitPad) this is redundant
/*
trex::THitSelection* trex::TTPCVolGroup::GetHitSelection(){
  trex::THitSelection* allHits = new trex::THitSelection();

  for(std::map<long, trex::TTPCUnitVolume*>::iterator el = fHitMap.begin(); el != fHitMap.end(); ++el){
    std::vector< trex::THandle<trex::TTPCHitPad> > hits = el->second->GetHits();
    for(std::vector< trex::THandle<trex::TTPCHitPad> >::iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt){
      trex::THandle<trex::TTPCHitPad> hit = *hitIt;

      allHits->AddHit(hit);
    };
  };

  return allHits;
  }*/

void trex::TTPCVolGroup::SetLean(int& leanVar, int leanVal, bool force){
  if(leanVar == leanVal) return;
  fIsClosed = false;
  // if conflict exists, revert to zero; otherwise overwrite
  if(!force && (leanVar*leanVal < 0)){
    leanVar = 0;
  }
  else{
    leanVar = leanVal;
  };
}

void trex::TTPCVolGroup::CheckHits(){
  // cells whos hits to check
  for(std::map<long, trex::TTPCUnitVolume*>::iterator volEl = fHitMap.begin(); volEl != fHitMap.end(); ++volEl){
    trex::TTPCUnitVolume* vol = volEl->second;
    if(!vol){
      std::cout << "WARNING:  empty trex::TTPCUnitVolume* found in TTPCVolGroup" << std::endl;
    }
    else{
      std::cout << "  Attempting to access trex::TTPCUnitVolume hits…" << std::endl;
      std::vector< trex::TTPCHitPad* > hits = vol->GetHits();
      std::cout << "  …success!" << std::endl;
      for(std::vector< trex::TTPCHitPad* >::iterator hitIt = vol->GetHitsBegin(); hitIt != vol->GetHitsEnd(); ++hitIt){
        trex::TTPCHitPad* nhit = *hitIt;
        if (!nhit){
          std::cout << "WARNING:  non trex::TTPCHitPad hit found in TTPCPathVolume" << std::endl;
        }
      }
    }
  }
}

//MDH
//Not used - just as well since these create new TTPCUnitVolumes which would need to be managed...
/*
void trex::TTPCVolGroup::AddPseudoHit(trex::TTPCUnitVolume* hit){
  fIsClosed = false;
  // add cell to map
  fHitMap[hit->GetID()] = hit;
}
void trex::TTPCVolGroup::AddNewPseudoHit(int x, int y, int z, float q){
  fIsClosed = false;
  // get cell's id
  long id = fLayout->Mash(x, y, z);
  std::map<long, trex::TTPCUnitVolume*>::iterator el = fHitMap.find(id);
  if (el == fHitMap.end()){
    // if it doesn't already exist, make a new cell at specified x, y and z id with specified charge
    trex::TTPCUnitVolume* cHit = new trex::TTPCUnitVolume();
    cHit->SetCell(x, y, z, -3, -3, -3, id);
    cHit->AddCharge(q);

    AddPseudoHit(cHit);
  }
  else{
    // if it does, just add specified charge to its charge
    el->second->AddCharge(q);
  };
}
void trex::TTPCVolGroup::AddNewPseudoGauss(int x, int y, int z, float q, int axis, float sigmaX, float sigmaY, float sigmaZ){
  fIsClosed = false;
  float sigmaI=1.;
  float sigmaJ=1.;
  // assign sigma two appropriate 2D dimensions depending on specified axis
  if(axis == 1){
    sigmaI = sigmaY;
    sigmaJ = sigmaZ;
  }
  else if(axis == 2){
    sigmaI = sigmaZ;
    sigmaJ = sigmaX;
  }
  else if(axis == 3){
    sigmaI = sigmaX;
    sigmaJ = sigmaY;
  };
  // ensure non-zero sigma
  sigmaI = std::max(.01f, sigmaI);
  sigmaJ = std::max(.01f, sigmaJ);

  int windowI = 1 + int(sigmaI);
  int windowJ = 1 + int(sigmaJ);

  const int sizeI = windowI*2 + 1;
  const int sizeJ = windowJ*2 + 1;

  // build unnormalised gaussian 
  float** base = new float* [sizeI];
  for(int i=0; i<sizeI; i++) base[i] = new float [sizeJ];

  float norm = 0;
  for(int i=-windowI; i<=windowI; i++){
    for(int j=-windowJ; j<=windowJ; j++){
      float val = exp(-.5*(float(i*i)/(sigmaI*sigmaI) + float(j*j)/(sigmaJ*sigmaJ)));
      base[i+windowI][j+windowJ] = val;
      norm += val;
    };
  };
  norm /= q;

  // add normalised gaussian 
  for(int i=-windowI; i<=windowI; i++){
    for(int j=-windowJ; j<=windowJ; j++){
      if(axis == 1) AddNewPseudoHit(x, y+i, z+j, base[i+windowI][j+windowJ]/norm);
      if(axis == 2) AddNewPseudoHit(x+j, y, z+i, base[i+windowI][j+windowJ]/norm);
      if(axis == 3) AddNewPseudoHit(x+i, y+j, z, base[i+windowI][j+windowJ]/norm);
    };
    delete[] base[i+windowI];
  };
  delete[] base;
}
*/
