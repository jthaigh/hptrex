#include "TTRExHVCluster.hxx"

void trex::TTRExHVCluster::CloseHits(){

  // Calculate the averages.
  fPosition.SetXYZ(0,0,0);
  fCharge = 0.0;

  // Calculate the charge weighted average position.
  double weightSum=0;
  for (auto h = fcHitPtrs.begin();
      h != fcHitPtrs.end();
      ++h) {

    trex::TTPCHitPad& hitPad = **h;

    double hitCharge = hitPad.GetCharge();
    double maxX = -99999.;
    double minX =  99999.;
    
    double tmpPadX = hitPad.GetPosition().X();
    if (tmpPadX > maxX) maxX = tmpPadX;
    if (tmpPadX < minX) minX = tmpPadX;
    fDeltaDrift = maxX-minX;

    weightSum += hitCharge;
    // warning Maybe use extrapolated Charge for the weighted position?
    fPosition[0] += hitCharge * hitPad.GetPosition().X();
    fPosition[1] += hitCharge * hitPad.GetPosition().Y();
    fPosition[2] += hitCharge * hitPad.GetPosition().Z();
    fCharge += hitCharge;

  }

  if (weightSum>0) {
    fPosition[0] = fPosition[0]/weightSum;
    fPosition[1] = fPosition[1]/weightSum;
    fPosition[2] = fPosition[2]/weightSum;
  }

}
