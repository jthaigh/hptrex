// eddy
#include "TTPCUnitVolume.hxx"

trex::TTPCUnitVolume::TTPCUnitVolume(){

  // initialise charge, hits and position to zero
  //fNegativePeakEarly = 0.;
  //fNegativePeakLate = 0.;
  fQ = 0.;
  fQMax = 0.;

  fHasPos = false;
  fFriendDist = 9999;

}

trex::TTPCUnitVolume::~TTPCUnitVolume(){
}

void trex::TTPCUnitVolume::SetCell(int x, int y, int z, long id){
  // set cell x, y and z id
  fX = x;
  fY = y;
  fZ = z;

  // set cell unique id
  fID = id;
}

void trex::TTPCUnitVolume::AddCharge(double q){
  // increment cell charge
  fQ += q;
  fQMax = std::max(fQMax, q);
}

void trex::TTPCUnitVolume::AddHit(trex::TTPCHitPad* hit){
  // add to list of hits
  fHits.push_back(hit);
}

void trex::TTPCUnitVolume::AddHits(std::vector< trex::TTPCHitPad* > hits){
  // associate multiple hits with this event
  for(std::vector< trex::TTPCHitPad* >::iterator hit = hits.begin(); hit != hits.end(); ++hit)
    AddHit(*hit);
}

void trex::TTPCUnitVolume::AddEvent(trex::TTPCHitPad* hit){
  // position, charge and time
  TVector3 pos = hit->GetPosition();
  double q = hit->GetCharge();

  // set charge weighted average position and time between old and input position
  fPos = ((fQ * fPos) + (q * pos)) * (1./(fQ + q));

  // increment charge
  AddCharge(q);

  // add the hit
  AddHit(hit);

  // indicate that cell now has a position
  fHasPos = true;
}


trex::TTPCCellInfo3D trex::TTPCUnitVolume::GetCellInfo3D(){
  trex::TTPCCellInfo3D padInfo;
  padInfo.x = fX;
  padInfo.y = fY;
  padInfo.z = fZ;

  return padInfo;
}
