#ifndef TTPCHitPad_hxx_seen
#define TTPCHitPad_hxx_seen

#include <vector>
#include "TVector3.h"
#include "TLorentzVector.h"

namespace trex {
  class TTPCHitPad;
}

/// Contains all the information relative to a waveform from a TPC pad.
/// This may contains more than one peak if the peaks are too closed
/// to be split into individual waveforms.
class trex::TTPCHitPad {
public:
  TTPCHitPad(double eDep, TLorentzVector pos4);

  double GetCharge(){return fChargeFit;}

  double GetTime(){return fTimeFit;}

  /// Short cut to access the Y position of the pad
  double Y(){return fPosition.Y();};
  /// Short cut to access the Z position of the pad
  double Z(){return fPosition.Z();};

  TVector3 GetPosition(){return fPosition;}

private:

  /// Amplitude for the peaks, as produced by the analytic gaussian fit (UNFIT if not fit)
  double fChargeFit;

  /// Fitted times for the peaks, as produced by the analytic gaussian fit (UNFIT if not fit)
  double fTimeFit;

  TVector3 fPosition;
};

#endif
