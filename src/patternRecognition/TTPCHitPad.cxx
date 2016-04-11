#include "TTPCHitPad.hxx"
#include <TGraph.h>
#include <TF1.h>


ClassImp(trex::TTPCHitPad);

//typedef std::vector< trex::THandle< trex::TSingleHit > >::const_iterator timebin;

trex::TTPCHitPad::TTPCHitPad(){
  InitParameters();
}


void trex::TTPCHitPad::InitParameters(){

}

void trex::TTPCHitPad::Init(){

/*  fNegativePeakBin = std::vector< trex::THandle<trex::TSingleHit> >();

  // Saturation
  fSaturation = 0;
  if ( this->GetTimeUncertainty() > 1.1* fSamplingTime) {
    fSaturation = int(this->GetTimeUncertainty() *2. / fSamplingTime);
  }

  fAtHoriEdge = trex::TGeomInfo::Get().TPC().PadIsAtHorEdge(this->GetGeomId());
  fAtVertEdge = trex::TGeomInfo::Get().TPC().PadIsAtVerEdge(this->GetGeomId());

  fHairCandidate = false;

  fDriftSense = int(trex::TGeomInfo::Get().TPC().GetDriftSense( this->GetGeomId() ));

  fChargeIntegral = 0.;
  trex::THandle<trex::TMultiHit> mh = GetMultiHit();
  int binNb = 0;
  double maxCharge = -9999.;
  for( timebin bin = mh->begin(); bin < mh->end(); bin++, binNb++){
    fChargeIntegral += (*bin)->GetCharge();
    if ((*bin)->GetCharge() > maxCharge){
      fMaxPeak = binNb;
      maxCharge = (*bin)->GetCharge();
    }
  }

  fIsFitted = kWFPEAKNOFIT;
  fChargeFit = -0xABCDEF;
  fTimeFit = -0xABCDEF;
  fSigmaFit = -0xABCDEF;

  // By default do the analytic fit.
  // The ROOT fit is available for studies but currently not considered
  // usable because probably too slow.
  AnalyticFit();*/
}
/*
unsigned int trex::TTPCHitPad::GetNumberPeaks(){
  return fPeakBin.size();
}


trex::THandle<trex::THit> trex::TTPCHitPad::GetMultiHit(){

  if ( this->GetContributorCount() > 0 ){
    return this->GetContributor(0);
  } else {
    return trex::THandle<trex::TMultiHit>();
  }
}


void trex::TTPCHitPad::AnalyticFit() {

  fIsFitted = kWFPEAKNOFIT;
  fChargeFit = -0xABCDEF;
  fTimeFit = -0xABCDEF;
  fSigmaFit = -0xABCDEF;

  // Only fit if there is only one peak in the hit.
  if (fPeakBin.size() > 1){
    return;
  }

  // Only perform analytic fit to peak if the peak is not saturated by more than 2 bins
  // Saturation at or below 2 bins is not distinguishable.
  if (fSaturation > 2){
    return;
  }

  // Only do the fit if there is enough room around the peak.
  trex::THandle<trex::TMultiHit> mh = GetMultiHit();
  if ((unsigned int)(fMaxPeak + fAnalyticFitRange) > mh->size() || (fMaxPeak - fAnalyticFitRange) < 0) {
    return;
  }

  // Check peak amplitude is large enough such that the fit can be
  // performed confidently to produce an accurate arrival time
  if (fCharge < fAnalyticFitAmpCut) {
    return;
  }

  // Primary electrons:
  double pe = 0.0;

  // Sum variables used in calculation of mean:
  double y = 0.0, // Sum_i y_i
         x = 0.0,      // Sum_i x_i
         x2 = 0.0,     // Sum_i x_i^2
         x3 = 0.0,     // Sum_i x_i^3
         x4 = 0.0,     // Sum_i x_i^4
         yx = 0.0,     // Sum_i y_i*x_i
         yx2 = 0.0,    // Sum_i y_i*x_i^2
         n = 0.0,      // Sum 1
         time,         // time bin
         time2;        // time bin squared

  typedef std::vector< trex::THandle< trex::TSingleHit > >::const_iterator iterator;
  iterator maxBin = mh->begin() + fMaxPeak;

  iterator start = maxBin - fAnalyticFitRange,
           end = maxBin + fAnalyticFitRange;

  // Calculation of sums:
  for (iterator j = start; j != end; j++) {
    pe = log((*j)->GetCharge());
    time = (*j)->GetTime();
    time2 = time*time;

    n += 1.0;
    y += pe;
    x += time;
    x2 += time2;
    x3 += time2*time;
    x4 += time2*time2;
    yx += pe*time;
    yx2 += pe*time2;
  }

  // Calculation of quantities obtained from matrix row reduction:
  double b = x2 - x*x/n,
         d = x3 - x2*x/n,
         c = yx - x*y/n,
         e = x3 - x2*x/n,
         f = x4 - x2*x2/n,
         g = yx2 - y*x2/n;

  // Calculation of mean, sigma, and amplitude based on a0, a1, and a2:
  double a2 = (g - e*c/b)/(f - e*d/b);
  double a1 = (c - a2*d)/b;
  double a0 = (y - a2*x2 - a1*x)/n;

  double analyticArrivalTime = -a1/(2*a2);
  double gaussianStDev = sqrt(-0.5/a2);
  double gaussianAmp = exp(a0 + analyticArrivalTime*a1 + analyticArrivalTime*analyticArrivalTime*a2);

  // Check success of fit based on the sigma of the gaussian. Inspection of fit
  // results indicates that when the fit fails a good indication of the failure
  // is a very large sigma, in excess of 1000.0 ns. Note the fitted location may
  // still be within a couple time bins. AG 8/7/2012
  if (gaussianStDev < 500.0) {
    fIsFitted = kWFPEAKANALYTICFIT;
    fChargeFit = gaussianAmp;
    fTimeFit = analyticArrivalTime;
    fSigmaFit = gaussianStDev;
  }
}


void trex::TTPCHitPad::RootGaussFit() {

  fIsFitted = kWFPEAKNOFIT;
  fChargeFit = -0xABCDEF;
  fTimeFit = -0xABCDEF;
  fSigmaFit = -0xABCDEF;

  if (fPeakBin.size() > 1){
    return;
  }

  std::vector<double> at;
  std::vector<double> aq;

  trex::THandle<trex::TMultiHit> mh = GetMultiHit();

  double localmaxADC = 0.0;
  int idxmax=-1;
  int idx=0;
  double   Time=0;

  typedef std::vector< trex::THandle< trex::TSingleHit > >::const_iterator iterator;
  for( iterator mhit = mh->begin(); mhit !=  mh->end(); mhit++ ) {
    if( localmaxADC < (*mhit)->GetCharge() ) {
      localmaxADC = (*mhit)->GetCharge();
      Time = (*mhit)->GetTime();
      idxmax = idx;
    }

    at.push_back( (*mhit)->GetTime() );
    aq.push_back( (*mhit)->GetCharge() );
    idx++;
  }

  // see if we can get better time estimate, and pulse width estimate
  int neachside=7;
  int nfit=std::min<int>(2*neachside+1,at.size());
  int idxmin = std::max<int>(0,idxmax-neachside);
  int idxmaxmax = std::min<int>(at.size()-1,idxmin+nfit);
  nfit = idxmaxmax-idxmin;

  TGraph tg(nfit, &(at[idxmin]), &(aq[idxmin]) );
  TF1 tf("tf","gaus",at[idxmin]-0.1,at[idxmaxmax]+0.1);
  tf.SetParameters( localmaxADC, Time, 100.0 ); // set a reasonable first guess
  tg.Fit(&tf,"");

  double Tfit = tf.GetParameter(1);
  double TSig = tf.GetParameter(2);
  double Qfit = tf.GetParameter(0)/0.133; //fudge factor to make it similar to integral charge

  // sanity check of fit.. if result yeilds
  // insane value, maybe reject this point?
  if ( ! (Tfit<at[idxmin] || Tfit>at[idxmaxmax] || Qfit<0.0) ){
    fIsFitted = kWFPEAKROOTGAUSSFIT;
    fChargeFit = Qfit;
    fTimeFit = Tfit;
    fSigmaFit = TSig;
  }

}
*/

std::vector<double> trex::TTPCHitPad::GetPeakTimes() {

  //std::vector<double> PeakTimes;
  //for (std::vector< trex::THandle< trex::TSingleHit > >::const_iterator bin = fPeakBin.begin(); bin !=  fPeakBin.end(); bin++) {
  //  PeakTimes.push_back((*bin)->GetTime());
  //}
  //return PeakTimes;
  std::vector<double> peakTimes;
  peakTimes.push_back(fTimeFit);
  return peakTimes;
}


std::vector<double> trex::TTPCHitPad::GetPeakCharges() {

  std::vector<double> PeakCharges;
//  for (std::vector< trex::THandle< trex::TSingleHit > >::const_iterator bin = fPeakBin.begin(); bin !=  fPeakBin.end(); bin++) {
//    PeakCharges.push_back((*bin)->GetCharge());
//  }
  PeakCharges.push_back(fChargeFit);
  return PeakCharges;
 
}
/*

std::vector<double> trex::TTPCHitPad::GetNegativePeakTimes() {

  std::vector<double> negativePeakTimes;
  for (std::vector< trex::THandle< trex::TSingleHit > >::const_iterator bin = fNegativePeakBin.begin(); bin !=  fNegativePeakBin.end(); bin++) {
    negativePeakTimes.push_back((*bin)->GetTime());
  }
  return negativePeakTimes;
}


std::vector<double> trex::TTPCHitPad::GetNegativePeakCharges() {

  std::vector<double> negativePeakCharges;
  for (std::vector< trex::THandle< trex::TSingleHit > >::const_iterator bin = fNegativePeakBin.begin(); bin !=  fNegativePeakBin.end(); bin++) {
    negativePeakCharges.push_back((*bin)->GetCharge());
  }
  return negativePeakCharges;
}


// void GetWaveformTimeRange(trex::THandle<trex::TReconHit> wf, double &min, double &max) {

//   // Simple single peak
//   if ( fPeakBin.size()){
//       min = fTime - fLowerTimeSpread;
//       max = fTime + fUpperTimeSpread;
//   // Many close peaks. Use the TSingleHits saved for each peak
//   } else {
//     min = wf->GetContributor(2)->GetTime() - fLowerTimeSpread;
//     max = wf->GetContributor(wf->GetContributorCount()-1)->GetTime() + fUpperTimeSpread;
//   }
// }

trex::THandle<trex::TMultiHit> trex::TTPCHitPad::ConvertToOAEvent() {
  return GetMultiHit();
}
*/
