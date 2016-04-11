// eddy
#include "TTPCUnitVolume.hxx"

const int trex::TTPCUnitVolume::scNTPC = 3;
const int trex::TTPCUnitVolume::scNHalf = 2;
const int trex::TTPCUnitVolume::scNMM = 12;
const int trex::TTPCUnitVolume::scNFEC = 6;
const int trex::TTPCUnitVolume::scNASIC = 4;

trex::TTPCUnitVolume::TTPCUnitVolume(){
  fID = 0;
  fTPC = 0;
  fHalf = 0;
  fMM = 0;
  fFEC = 0;
  fASIC = 0;
  fASICRegionY = -1;
  fASICRegionZ = -1;
  fUniqueASICID = -1;

  // initialise charge, hits and position to zero
  fSegX = fSegY = fSegZ = 0;
  fNegativePeakEarly = 0.;
  fNegativePeakLate = 0.;
  fQ = 0.;
  fQMax = 0.;

  fDeltaTagged = false;
  fLowChargeTagged = false;
  fEarlyNegativeTagged = false;
  fLateNegativeTagged = false;
  fFullASICTagged = false;
  fSatASICTagged = false;

  fTime = 0.;
  fTimeNom = 0.;
  fTimeMin = +99999999.;
  fTimeMax = -99999999.;
  fHasPos = false;
  fFriendDist = 9999;

}

trex::TTPCUnitVolume::~TTPCUnitVolume(){
}

void trex::TTPCUnitVolume::SetCell(int x, int y, int z, int edgeX, int edgeY, int edgeZ, long id){
  // set cell x, y and z id
  fX = x;
  fY = y;
  fZ = z;

  // set cell x, y and z edge status
  fEdgeX = edgeX;
  fEdgeY = edgeY;
  fEdgeZ = edgeZ;

  // set cell unique id
  fID = id;
}

void trex::TTPCUnitVolume::SetAux(int segX, int segY, int segZ){
  // set cell x, y and z id
  fSegX = segX;
  fSegY = segY;
  fSegZ = segZ;
}

void trex::TTPCUnitVolume::SetRegion(int asicRegionY, int asicRegionZ){
  fASICRegionY = asicRegionY;
  fASICRegionZ = asicRegionZ;
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
  double time = 0.;
  double timeNom = 0.;

  // set average, max and min times from hit peak times
  std::vector<double> peakTimes = hit->GetPeakTimes();
  if(peakTimes.size()){
    time = 0.;
    for(std::vector<double>::iterator peakTimeIt = peakTimes.begin(); peakTimeIt != peakTimes.end(); ++peakTimeIt){
      double offPeakTime = *peakTimeIt - fTimeOffset;
      fTimeMin = std::min(fTimeMin, offPeakTime);
      fTimeMax = std::max(fTimeMax, offPeakTime);
      timeNom += offPeakTime;
      time += *peakTimeIt;
    };
    timeNom /= (double)peakTimes.size();
    time /= (double)peakTimes.size();
  }
  else{
    double offPeakTime = hit->GetTime() - fTimeOffset;
    fTimeMin = std::min(fTimeMin, offPeakTime);
    fTimeMax = std::max(fTimeMax, offPeakTime);
    timeNom = offPeakTime;
    time = hit->GetTime();
  };

  // set charge weighted average position and time between old and input position
  fPos = ((fQ * fPos) + (q * pos)) * (1./(fQ + q));
  fTimeNom = ((fQ * fTimeNom) + (q * timeNom)) * (1./(fQ + q));
  fTime = ((fQ * fTime) + (q * time)) * (1./(fQ + q));

  // also add negative peak charges before and after
  std::vector<double> negativePeakCharges = hit->GetNegativePeakCharges();
  std::vector<double> negativePeakTimes = hit->GetNegativePeakTimes();
  unsigned int negativePeakN = negativePeakCharges.size();
  for(unsigned int i=0; i<negativePeakN; ++i){
    double negativeCharge = negativePeakCharges.at(i);
    double negativeTime = negativePeakTimes.at(i);

    if(negativeTime < time){
      fNegativePeakEarly = std::min(fNegativePeakEarly, negativeCharge);
    }
    else if(negativeTime > time){
      fNegativePeakLate = std::min(fNegativePeakLate, negativeCharge);
    };
  };

  // increment charge
  AddCharge(q);

  // add the hit
  AddHit(hit);

  // indicate that cell now has a position
  fHasPos = true;
}

void trex::TTPCUnitVolume::GetEdges(int& edgeX, int& edgeY, int& edgeZ){
  // return cell x, y and z edge status
  edgeX = fEdgeX;
  edgeY = fEdgeY;
  edgeZ = fEdgeZ;
}

trex::TTPCCellInfo3D trex::TTPCUnitVolume::GetCellInfo3D(){
  trex::TTPCCellInfo3D padInfo;
  padInfo.x = fX;
  padInfo.y = fY;
  padInfo.z = fZ;
  padInfo.edgeX = fX;
  padInfo.edgeY = fY;
  padInfo.edgeZ = fZ;

  return padInfo;
}

int trex::TTPCUnitVolume::GetNPeaksSum(){
  unsigned int peaksSum=0;
  for(std::vector< trex::TTPCHitPad* >::iterator hitIt = fHits.begin(); hitIt != fHits.end(); ++hitIt){
    trex::TTPCHitPad* hitPad = *hitIt;
    peaksSum += hitPad->GetNumberPeaks();
  }
  return (int)peaksSum;
}
int trex::TTPCUnitVolume::GetNPeaksMax(){
  unsigned int peaksSum=0;
  for(std::vector< trex::TTPCHitPad* >::iterator hitIt = fHits.begin(); hitIt != fHits.end(); ++hitIt){
    trex::TTPCHitPad* hitPad = *hitIt;
    peaksSum = std::max(peaksSum, hitPad->GetNumberPeaks());
  }
  return (int)peaksSum;
}
int trex::TTPCUnitVolume::GetNSaturated(){
  int nSat=0;
  for(std::vector< trex::TTPCHitPad* >::iterator hitIt = fHits.begin(); hitIt != fHits.end(); ++hitIt){
    trex::TTPCHitPad* hitPad = *hitIt;
    nSat += (int)(hitPad->Saturation() > 1);
  };
  return nSat;
}
int trex::TTPCUnitVolume::GetSaturation(){
  int satSum=0;
  for(std::vector< trex::TTPCHitPad* >::iterator hitIt = fHits.begin(); hitIt != fHits.end(); ++hitIt){
    trex::TTPCHitPad* hitPad = *hitIt;
    satSum += hitPad->Saturation();
  };
  return satSum;
}

void trex::TTPCUnitVolume::DefineUniqueASICID(){
  long tpcPart = scNASIC * scNFEC * scNMM * scNHalf * fTPC;
  long halfPart = scNASIC * scNFEC * scNMM * fHalf;
  long mmPart = scNASIC * scNFEC * fMM;
  long fecPart = scNASIC * fFEC;
  long asicPart = fASIC;

  fUniqueASICID = tpcPart + halfPart + mmPart + fecPart + asicPart;
}
