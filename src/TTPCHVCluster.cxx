#include <cmath>
#include <algorithm>

#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoNode.h>

#include <HEPUnits.hxx>
#include <THitSelection.hxx>
#include <TOADatabase.hxx>
#include <TGeomIdManager.hxx>
#include <TRecPackManager.hxx>
#include <TOARuntimeParameters.hxx>

#include "TTPCHVCluster.hxx"
#include "TTPCHitPad.hxx"
#include "TTPCCalibration.hxx"

// *********************************************************************************
ND::TTPCHVCluster::TTPCHVCluster()
  : fHitsOpen(true), fGeomId(0), fCharge(0.0), fIntCharge(0.0), fTime(0.0),
      fDeltaDrift(0.0), fPosition(0,0,0), fDriftDistance(0.0), fT0(0.0),
      fSpread(100*unit::meter,100*unit::meter,100*unit::meter),
      fUncertainty(100*unit::meter,100*unit::meter,100*unit::meter),
      fTimeUncertainty(1*unit::ns),
      fIsXHit(false), fIsYHit(false), fIsZHit(false),
      fIsVertical(true), fEndNode(false), fFieldCorrection(0.0),
      fMaxCharge(0.0), fMaxNPeaks(0), fNSaturated(0), 
      fAtVertEdge(false), fAtHoriEdge(false), fSuspectPadTiming(false),
      fAllFitted(true), fOkForSeed(true), fOkForFit(true), fUsable(true){
  fHits.SetTitle("TPC Cluster Hits");
  SetBit(kCanDelete,true);
  fDeltaYZ[0] = 0.0;
  fDeltaYZ[1] = 0.0;

  fUseExtrapolatedCharge = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.Cluster.UseExtrapolatedCharge");                       // Use extrapolated charge for saturated hits?
  fUseFittedCharge = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.Cluster.UseFittedCharge");                                   // Use fitted charge when present?
  fUseExtrapolatedIntCharge = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.Cluster.UseExtrapolatedIntegratedCharge");          // Use extrapolated integrated charge for saturated hits?
  fHeightToAreaFactor = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.Cluster.ExtrapolatedHeightToIntegralConversionFactor");   // Factor to convert extrapolated heights to extrapolated integrals
}
// *********************************************************************************
ND::TTPCHVCluster::~TTPCHVCluster() { }


// *********************************************************************************
ND::TTPCHVCluster::TTPCHVCluster(ND::THitSelection& hits, bool isAVerticalCluster)
  : fHitsOpen(true), fGeomId(0), fCharge(0.0), fIntCharge(0.0), fTime(0.0),
      fDeltaDrift(0.0), fPosition(0,0,0), fDriftDistance(0.0), fT0(0.0),
      fSpread(100*unit::meter,100*unit::meter,100*unit::meter),
      fUncertainty(100*unit::meter,100*unit::meter,100*unit::meter),
      fTimeUncertainty(1*unit::ns),
      fIsXHit(false), fIsYHit(false), fIsZHit(false),
      fIsVertical(isAVerticalCluster), fEndNode(false), fFieldCorrection(0.0),
      fMaxCharge(0.0), fMaxNPeaks(0), fNSaturated(0), 
      fAtVertEdge(false), fAtHoriEdge(false), fSuspectPadTiming(false),
      fAllFitted(true), fOkForSeed(true), fOkForFit(true), fUsable(true){
  fHits.SetTitle("TPC Cluster Hits");
  SetBit(kCanDelete,true);
  fDeltaYZ[0] = 0.0;
  fDeltaYZ[1] = 0.0;

  fUseExtrapolatedCharge = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.Cluster.UseExtrapolatedCharge");                       // Use extrapolated charge for saturated hits?
  fUseFittedCharge = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.Cluster.UseFittedCharge");                                   // Use fitted charge when present?
  fUseExtrapolatedIntCharge = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.Cluster.UseExtrapolatedIntegratedCharge");          // Use extrapolated integrated charge for saturated hits?
  fHeightToAreaFactor = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.Cluster.ExtrapolatedHeightToIntegralConversionFactor");   // Factor to convert extrapolated heights to extrapolated integrals

  AddHitSelection(hits);
  CloseHits();
}
  


// *********************************************************************************
ND::TGeometryId ND::TTPCHVCluster::GetGeomId(void) const {
    if (fHitsOpen) CloseHits();
    return fGeomId;
}

// *********************************************************************************
const TVector3& ND::TTPCHVCluster::GetPosition(void) const {
    if (fHitsOpen) CloseHits();
    return fPosition;
}

// *********************************************************************************
const TVector3& ND::TTPCHVCluster::GetCalibPosition(void) const {
    if (fHitsOpen) CloseHits();
    return fCalibratedPosition;
}

// *********************************************************************************
double ND::TTPCHVCluster::GetCharge(void) const {
    if (fHitsOpen) CloseHits();
    return fCharge;
}

// *********************************************************************************
double ND::TTPCHVCluster::GetMaxCharge(void) const {
    if (fHitsOpen) CloseHits();
    return fMaxCharge;
}

// *********************************************************************************
double ND::TTPCHVCluster::GetIntCharge(void) const {
  if (fHitsOpen) CloseHits();
  return fIntCharge;
}

// *********************************************************************************
double ND::TTPCHVCluster::GetTime(void) const {
    if (fHitsOpen) CloseHits();
    return fTime;
}

// *********************************************************************************
double ND::TTPCHVCluster::GetDeltaDrift(void) const {
    if (fHitsOpen) CloseHits();
    return fDeltaDrift;
}


// *********************************************************************************
const TVector3& ND::TTPCHVCluster::GetSpread(void) const {
    if (fHitsOpen) CloseHits();
    return fSpread;
}

// *********************************************************************************
const TVector3& ND::TTPCHVCluster::GetUncertainty(void) const {
    if (fHitsOpen) CloseHits();
    return fUncertainty;
}

// *********************************************************************************
double ND::TTPCHVCluster::GetTimeUncertainty(void) const {
    if (fHitsOpen) CloseHits();
    return fTimeUncertainty;
}

// *********************************************************************************
bool ND::TTPCHVCluster::IsXHit(void) const {
    if (fHitsOpen) CloseHits();
    return fIsXHit;
}

// *********************************************************************************
bool ND::TTPCHVCluster::IsYHit(void) const {
    if (fHitsOpen) CloseHits();
    return fIsYHit;
}

// *********************************************************************************
bool ND::TTPCHVCluster::IsZHit(void) const {
    if (fHitsOpen) CloseHits();
    return fIsZHit;
}

// *********************************************************************************
TGeoNode* ND::TTPCHVCluster::GetParentNode(int i) const {
    if (fHitsOpen) CloseHits();
    if (i<0) return NULL;
    if (i>=(int) fParentNodes.size()) return fParentNodes.back();
    return fParentNodes[i];
}

// *********************************************************************************
const ND::THitSelection& ND::TTPCHVCluster::GetHits(void) const {
    return fHits;
}

//////////////////////////////////////////////////
// Setter methods to add hits to the TComboHit.
//////////////////////////////////////////////////
// *********************************************************************************
void ND::TTPCHVCluster::AddHit(ND::THandle<ND::THit>& hit) {
    OpenHits();
    fHits.push_back(hit);
}

// *********************************************************************************
void ND::TTPCHVCluster::AddHitSelection(ND::THitSelection& hits) {
    for (ND::THitSelection::iterator h = hits.begin();
         h != hits.end();
         ++h) {
        AddHit(*h);
    }
}

// *********************************************************************************
void ND::TTPCHVCluster::OpenHits() {
    fHitsOpen = true;
}

// *********************************************************************************
void ND::TTPCHVCluster::CloseHits() const {
  if (!fHitsOpen) return;
  if (fHits.size()<1) return;
  fHitsOpen = false;

  // Calculate the averages.
  fPosition.SetXYZ(0,0,0);
  fCharge = 0.0;
  fIntCharge = 0.0;
  fTime = 0.0;
  fSpread.SetXYZ(1,1,1);
  fUncertainty.SetXYZ(0,0,0);
  fTimeUncertainty = 1.0*unit::ns;
  fIsXHit = false;
  fIsYHit = false;
  fIsZHit = false;


  // Calculate the charge weighted average position.
  double weightSum=0;
  for (ND::THitSelection::const_iterator h = fHits.begin();
      h != fHits.end();
      ++h) {
    ND::THandle<ND::TTPCHitPad> hitPad = (*h);
    if (!hitPad){
      std::cout<<"TTPCHVCluster ERROR: A hit in the cluster is not a TTPCHitPad !"<<std::endl;
      continue;
    }

    double hitCharge = hitPad->GetCharge();
    double hitIntCharge = hitPad->ChargeIntegral();
    double hitTime = hitPad->GetTime();
    double maxX = -99999.;
    double minX =  99999.;
    if (fUseFittedCharge and hitPad->IsFitted()){ // Use fitted charge?
      hitCharge = hitPad->ChargeFit();
      hitTime = hitPad->TimeFit();
      if (fabs(hitPad->TimeFit() - hitPad->GetTime()) > 100)
        fSuspectPadTiming = true;
    } else {
      fAllFitted = false;
      if (fUseExtrapolatedCharge and hitPad->Saturation() > 2) { // Use extrapolated charge?
        hitCharge = hitPad->ChargeExtrapolated();
      }
    }
    if (fUseExtrapolatedIntCharge and hitPad->Saturation() > 2) { // Use extrapolated integrated charge?
      // Make sure the extrapolated charge is actually larger than the measured one
      double tempCharge = hitPad->ChargeExtrapolated() * fHeightToAreaFactor;
      if (tempCharge > hitIntCharge){
        hitIntCharge = tempCharge;
      }
    }

    fTime += hitCharge * hitTime;
    double tmpPadX = hitTime * ND::tpcCalibration().GetDriftVelocity();
    if (tmpPadX > maxX) maxX = tmpPadX;
    if (tmpPadX < minX) minX = tmpPadX;
    fDeltaDrift = maxX-minX;

    weightSum += hitCharge;
    // warning Maybe use extrapolated Charge for the weighted position?
    fPosition[0] += hitCharge * (*h)->GetPosition().X();
    fPosition[1] += hitCharge * (*h)->GetPosition().Y();
    fPosition[2] += hitCharge * (*h)->GetPosition().Z();
    fCharge += hitCharge;
    fIntCharge += hitIntCharge;

    if( hitCharge > fMaxCharge) fMaxCharge = hitCharge;
    fAtVertEdge = fAtVertEdge || hitPad->IsAtVertEdge();
    fAtHoriEdge = fAtHoriEdge || hitPad->IsAtHoriEdge();
    // Count the number of saturated waveforms.
    // 2 bins of same height can happen,
    // so only when 3 or more bins have the max ADC is the waveform considered saturated.
    if( hitPad->Saturation() > 2 ) fNSaturated++;
  }

  if (weightSum>0) {
    fPosition[0] = fPosition[0]/weightSum;
    fPosition[1] = fPosition[1]/weightSum;
    fPosition[2] = fPosition[2]/weightSum;
    fTime = fTime/weightSum;
  }

  // Calculate the "calibrated" position.
  // Position of the pad that measured the charge, i.e. the read out plane X.
  double RPx = this->GetPosition().X();
  ND::THandle<ND::TTPCHitPad> hitPad = *(fHits.begin());
  double Sense = hitPad->DriftSense();
  fDriftDistance = (fTime - fT0 - ND::tpcCalibration().GetTimeOffset()) * ND::tpcCalibration().GetDriftVelocity();
  // How far are we from the read out plane ?
  fCalibratedPosition.SetX( RPx - (Sense * fDriftDistance));
  fCalibratedPosition.SetY( this->GetPosition().Y() + fDeltaYZ[0]);
  fCalibratedPosition.SetZ( this->GetPosition().Z() + fDeltaYZ[1]);

  // =====================================================================================
  // From this line of = to the next, the code is just the copy from TComboHit::CloseHits.
  // It's not clear what we want for the uncertainties and the spread here.
  // The spread of the hits.
  double sprX=0, sprY=0, sprZ=0;
  // The minimum spread of any hit in the collection.
  double minSprX=1000*unit::mm, minSprY=1000*unit::mm, minSprZ=1000*unit::mm;
  // The position uncertainties.
  double uncXX=0, uncYY=0, uncZZ=0; 
  double uncXY=0, uncXZ=0, uncYZ=0; 
  // The weighting for the position uncertainties.
  double wghXX=0, wghYY=0, wghZZ=0; 
  double wghXY=0, wghXZ=0, wghYZ=0; 
  for (ND::THitSelection::const_iterator h = fHits.begin();
      h != fHits.end();
      ++h) {
    double dx = (*h)->GetPosition().X() - fPosition.X();
    double dy = (*h)->GetPosition().Y() - fPosition.Y();
    double dz = (*h)->GetPosition().Z() - fPosition.Z();
    double dt = (*h)->GetTime() - fTime;
    double sx = 1.0/sqrt(fCharge);
    double sy = 1.0/sqrt(fCharge);
    double sz = 1.0/sqrt(fCharge);

    uncXX += dx*dx/(sx*sx);
    wghXX += 1.0/(sx*sx);
    uncYY += dy*dy/(sy*sy);
    wghYY += 1.0/(sy*sy);
    uncZZ += dz*dz/(sz*sz);
    wghZZ += 1.0/(sz*sz);
    uncXY += dx*dy/(sx*sy);
    wghXY += 1.0/(sx*sy);
    uncXZ += dx*dz/(sx*sz);
    wghXZ += 1.0/(sx*sz);
    uncYZ += dy*dz/(sy*sz);
    wghYZ += 1.0/(sy*sz);
    fTimeUncertainty += (*h)->GetCharge()*dt*dt;
    sprX = std::max((*h)->GetSpread().X(),sprX);
    sprX = std::max(dx,sprX);
    sprY = std::max((*h)->GetSpread().Y(),sprY);
    sprY = std::max(dy,sprY);
    sprZ = std::max((*h)->GetSpread().Z(),sprZ);
    sprZ = std::max(dz,sprZ);
    minSprX = std::min(minSprX, (*h)->GetSpread().X());
    minSprY = std::min(minSprY, (*h)->GetSpread().Y());
    minSprZ = std::min(minSprZ, (*h)->GetSpread().Z());
  }
  fSpread.SetXYZ(sprX,sprY,sprZ);

  // The uncXX .. uncYZ variables will hold the average width of the
  // distribution [i.e. avg(x - avg(x))]
  if (fCharge>0) {
    uncXX /= wghXX;
    uncYY /= wghYY;
    uncZZ /= wghZZ;
    uncXY /= wghXY;
    uncXZ /= wghXZ;
    uncYZ /= wghYZ;
    fTimeUncertainty /= fCharge;
  }

  // Make sure that width of the distribution isn't smaller than the
  // smallest geometric element.  The "4" comes since minSpr is half of the
  // element size, and the "12" comes from assuming the real hit position is
  // spread uniformly over the element.
  uncXX = std::max(uncXX,4*minSprX*minSprX/12);
  uncYY = std::max(uncYY,4*minSprY*minSprY/12);
  uncZZ = std::max(uncZZ,4*minSprZ*minSprZ/12);

  // The uncXX .. uncYZ variables will now hold the variance of the average
  // position.  Normalize for the number of measurements where one unit of
  // charge is treated as one measurement.  That is the right thing for
  // photosensor based detectors (as long as the charge is measured in
  // photo-electrons), but I'm not sure if it is OK for TPC based
  // measurements.  This is separate from the weighting above which is used
  // to calculate the average width of the distribution.  The normalization
  // above accounts for the multiple measurements of the position.
  double norm = std::max(double(fCharge-1), double(fHits.size()-1));
  if (norm>0) {
    uncXX /= norm;
    uncYY /= norm;
    uncZZ /= norm;
    uncXY /= norm;
    uncXZ /= norm;
    uncYZ /= norm;
    fTimeUncertainty /= norm;
  }

  // Set the uncertainty to the width of the distribution.
  fUncertainty.SetXYZ(sqrt(uncXX),sqrt(uncYY),sqrt(uncZZ));
  fTimeUncertainty = sqrt(fTimeUncertainty);
  // =================================================================================

  // Use the first node as reference just because.
  TVector3 firstHitPos = (*fHits.begin())->GetPosition();
  ND::TOADatabase::Get().GeomId().GetGeometryId(
      firstHitPos.X(), firstHitPos.Y(), firstHitPos.Z(), fGeomId);


  // TODO: Modify to do things more TPC-like
  // For example for vertical clusters, use the charged weighted Y position for the cluster position.
  // The spread could also be defined to be useful for a chi2 calculation.
}



// *********************************************************************************
void ND::TTPCHVCluster::PreparePropagSurf(const std::string surfname){
  // Don't create a new surface if there is already one available.
  if (fPropagSurfName.length() > 0)
    return;
  // The propagation surface is defined by the cluster orientation.
  // Maybe for large X angle we should switch to YZ planes ...
  TVector3 projNorm(0.,0.,0.);
  if ( IsVertical() ){
    projNorm.SetZ(1.);
  } else {
    projNorm.SetY(1.);
  }
  // Build the surface
  bool ok = ND::tman().BuildPropagationSurface(fCalibratedPosition, projNorm, fPropagSurf);
  if (!ok){
    std::cout<<"TTPCHVCluster: ERROR: Couldn't build the propagation surface !"<<std::endl;
    //TODO: proper exception
    throw;
  }

  fPropagSurfName = surfname;
  const std::string volname = fPropagSurf.name(RP::setup_name);
  ok = ND::rpman().geometry_svc().setup().add_surface(volname,surfname, &fPropagSurf);
  if (!ok){
    std::cout<<"TTPCHVCluster: ERROR: Couldn't add the propagation surface to the RecPack geometry !"<<std::endl;
    //TODO: proper exception
    throw;
  }
}


// *********************************************************************************
Surface& ND::TTPCHVCluster::GetPropagSurf(){
  return fPropagSurf;
}


// *********************************************************************************
void ND::TTPCHVCluster::RemovePropagSurf(){
  ND::rpman().geometry_svc().setup().remove_surface(fPropagSurfName);
}


// *********************************************************************************
void ND::TTPCHVCluster::SetIsHorizontal(void){
  fIsVertical = false;
}


// *********************************************************************************
void ND::TTPCHVCluster::SetIsVertical(void){
  fIsVertical = true;
}


// *********************************************************************************
bool ND::TTPCHVCluster::IsHorizontal(void){
  return !fIsVertical;
}


// *********************************************************************************
bool ND::TTPCHVCluster::IsVertical(void){
  return fIsVertical;
}


// *********************************************************************************
void ND::TTPCHVCluster::SetEndNode(void){
  fEndNode = true;
}


// *********************************************************************************
void ND::TTPCHVCluster::ClearEndNode(void){
  fEndNode = false;
}


// *********************************************************************************
bool ND::TTPCHVCluster::IsEndNode(void){
  return fEndNode;
}


// *********************************************************************************
double ND::TTPCHVCluster::GetFieldCorrection(void){
  return fFieldCorrection;
}


// *********************************************************************************
double ND::TTPCHVCluster::GetFieldCorrection(void) const{
  return fFieldCorrection;
}


// *********************************************************************************
void ND::TTPCHVCluster::SetFieldCorrection(double val){
  fFieldCorrection = val;
}


// *********************************************************************************
bool ND::TTPCHVCluster::IsAtVertEdge(){
  return fAtVertEdge;
}


// *********************************************************************************
bool ND::TTPCHVCluster::IsAtHoriEdge(){
  return fAtHoriEdge;
}


// *********************************************************************************
bool ND::TTPCHVCluster::HasSuspiciousPadTiming(){
  return fSuspectPadTiming;
}


// *********************************************************************************
unsigned int ND::TTPCHVCluster::GetMaxNPeaks(){
  return fMaxNPeaks;
}


// *********************************************************************************
unsigned int ND::TTPCHVCluster::GetNSaturated(){
  return fNSaturated;
}


// *********************************************************************************
int ND::TTPCHVCluster::AllFitted(){
  return fAllFitted;
}


// *********************************************************************************
void ND::TTPCHVCluster::SetT0(double T0){
  fT0 = T0;
  // position of the pad that measured the charge, i.e. the read out plane X.
  double RPx = this->GetPosition().X();
  ND::THandle<ND::TTPCHitPad> hitPad = *(fHits.begin());
  double Sense = hitPad->DriftSense();
  // How far are we from the read out plane ?
  fDriftDistance = (fTime - fT0 - ND::tpcCalibration().GetTimeOffset()) * ND::tpcCalibration().GetDriftVelocity();
  fCalibratedPosition.SetX( RPx - (Sense * fDriftDistance));
}

// *********************************************************************************
void ND::TTPCHVCluster::SetDeltaYZ(double dY, double dZ){
  fDeltaYZ[0] = dY;
  fDeltaYZ[1] = dZ;
}

// *********************************************************************************
double ND::TTPCHVCluster::GetDeltaY(){
  return fDeltaYZ[0];
}

// *********************************************************************************
double ND::TTPCHVCluster::GetDeltaZ(){
  return fDeltaYZ[1];
}

// *********************************************************************************
double ND::TTPCHVCluster::X(){
  return this->GetPosition().X();
}

// *********************************************************************************
double ND::TTPCHVCluster::Y(){
  return this->GetPosition().Y();
}

// *********************************************************************************
double ND::TTPCHVCluster::Z(){
  return this->GetPosition().Z();
}

// *********************************************************************************
TVector3 ND::TTPCHVCluster::GetCalibratedPosition(){
  return fCalibratedPosition;
}

// *********************************************************************************
double ND::TTPCHVCluster::GetDriftDistance(){
  return fDriftDistance;
}

// *********************************************************************************
double ND::TTPCHVCluster::CalibX(){
  return fCalibratedPosition.X();
}

// *********************************************************************************
double ND::TTPCHVCluster::CalibY(){
  return fCalibratedPosition.Y();
}

// *********************************************************************************
double ND::TTPCHVCluster::CalibZ(){
  return fCalibratedPosition.Z();
}

// *********************************************************************************
bool ND::TTPCHVCluster::hasT0(){
  return fT0;
}

// *********************************************************************************
bool ND::TTPCHVCluster::hasDeltaYZ(){
  return (fDeltaYZ[0] || fDeltaYZ[1]);
}


// *********************************************************************************
ND::THandle<ND::TComboHit> ND::TTPCHVCluster::ConvertToOAEvent() {
  ND::THandle<ND::TComboHit> combo( new ND::TComboHit ); 
  for (ND::THitSelection::const_iterator h = fHits.begin(); h != fHits.end(); ++h) {
    ND::THandle<ND::TTPCHitPad> tmpHP = (*h);
    ND::THandle<ND::TMultiHit> tmpMH = tmpHP->ConvertToOAEvent();
    ND::THandle<ND::THit> tmpHit = tmpMH;
    combo->AddHit(tmpHit);
  }
  return combo;
}


