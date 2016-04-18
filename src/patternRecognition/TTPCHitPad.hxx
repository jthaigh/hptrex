#ifndef TTPCHitPad_hxx_seen
#define TTPCHitPad_hxx_seen

#include <vector>
#include "TVector3.h"
#include "TLorentzVector.h"

/// Possible sources of the T0
enum {kWFPEAKNOFIT = 0, kWFPEAKANALYTICFIT, kWFPEAKROOTGAUSSFIT};

namespace trex {
  class TTPCHitPad;
}

/// Contains all the information relative to a waveform from a TPC pad.
/// This may contains more than one peak if the peaks are too closed
/// to be split into individual waveforms.
class trex::TTPCHitPad {
public:
  TTPCHitPad(double eDep, TLorentzVector pos4);
  /// Returns a vector with the charge of each one the peaks in the waveform.
  /// If there is only one peak, the given charge is the same as GetCharge().
  std::vector<double> GetPeakCharges();
  double GetCharge(){return GetPeakCharges()[0];}

  /// Returns a vector with the time of each one the peaks in the waveform.
  /// If there is only one peak, the given time is the same as GetTime().
  std::vector<double> GetPeakTimes();
  double GetTime(){return GetPeakTimes()[0];}

  std::vector<double> GetNegativePeakTimes(){std::vector<double> dummy;return dummy;}
  std::vector<double> GetNegativePeakCharges(){std::vector<double> dummy;return dummy;}

  unsigned int GetNumberPeaks(){return 1;}
  /// Return a vector with the charges of each of the negative peaks in the waveform
  //  std::vector<double> GetNegativePeakCharges();
  /// Return a vector with the times of each of the negative peaks in the waveform
  //std::vector<double> GetNegativePeakTimes();

  /// Number of bins with the same value at the highest peak.
  /// A return value of 1 is for a peak considered not saturated.
  int Saturation(){return fSaturation;};
  int DriftSense(){return fDriftSense;};

  /// Returns the charge integral for the waveform
  /// as opposed to the height of the highest peak returned by GetCharge()
  double ChargeIntegral(){return fChargeIntegral;};

  /// Returns the extrapolated charge of the highest peak in the waveform
  /// as opposed to the possibly saturated height returned by GetCharge()
  double ChargeExtrapolated(){return fChargeExtrapolated;};

  // Is this pad against the vertical edge of the MM ?
  bool IsAtVertEdge(){return fAtVertEdge;};
  // Is this pad against the horizontal edge of the MM ?
  bool IsAtHoriEdge(){return fAtHoriEdge;};

  /// Report if the peak of the waveform has been fitted and by which algorithm:
  /// kWFPEAKNOFIT        -> No successful fit
  /// kWFPEAKANALYTICFIT  -> Analytic fit
  /// kWFPEAKROOTGAUSSFIT -> ROOT Gaussian fit
  int IsFitted(){return fIsFitted;};

  /// Returns the charge obtained from the fit to the peak
  double ChargeFit(){return fChargeFit;};
  /// Returns the time obtained from the fit to the peak
  double TimeFit(){return fTimeFit;};
  /// Returns the sigma obtained from the fit to the peak
  double SigmaFit(){return fSigmaFit;};

  /// Short cut to access the Y position of the pad
  double Y(){return fPosition.Y();};
  /// Short cut to access the Z position of the pad
  double Z(){return fPosition.Z();};

  TVector3 GetPosition(){return fPosition;}

private:
  /// Initialise member variables
  void Init();
  void InitParameters();

  /// Index of the bin of maximum charge in the waveform
  int fMaxPeak;
  unsigned int fSaturation;
  int fDriftSense;

  bool fAtHoriEdge;
  bool fAtVertEdge;

  double fChargeIntegral;
  double fChargeExtrapolated;

  int fIsFitted;

  /// Amplitude for the peaks, as produced by the analytic gaussian fit (UNFIT if not fit)
  double fChargeFit;
  /// Fitted times for the peaks, as produced by the analytic gaussian fit (UNFIT if not fit)
  double fTimeFit;
  /// Gaussian sigma from the peak fit (UNFIT if not fit)
  double fSigmaFit;

  /// tpcRecon.Reco.Wave.AnalyticFitRange
  int fAnalyticFitRange;
  /// tpcRecon.Reco.Wave.AnalyticFitAmpCut
  double fAnalyticFitAmpCut;
  /// tpcRecon.Reco.Wave.LowerTimeSpread
  double fLowerTimeSpread;
  /// tpcRecon.Reco.Wave.UpperTimeSpread
  double fUpperTimeSpread;
  /// tpcRecon.AfterTPC.SamplingTime
  double fSamplingTime;

  TVector3 fPosition;
};

#endif
