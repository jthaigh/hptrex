#ifndef TTPCHitPad_hxx_seen
#define TTPCHitPad_hxx_seen

#include <TReconHit.hxx>
#include <TMultiHit.hxx>

/// Possible sources of the T0
enum {kWFPEAKNOFIT = 0, kWFPEAKANALYTICFIT, kWFPEAKROOTGAUSSFIT};

namespace ND {
  class TTPCHitPad;
}

/// Contains all the information relative to a waveform from a TPC pad.
/// This may contains more than one peak if the peaks are too closed
/// to be split into individual waveforms.
class ND::TTPCHitPad : public TReconHit {
public:
  TTPCHitPad();
  TTPCHitPad(const TWritableReconHit& val, std::vector< ND::THandle< ND::TSingleHit > >, double ChargeExtrapolated);
  TTPCHitPad(const TWritableReconHit& val, std::vector< ND::THandle< ND::TSingleHit > > bins, std::vector< ND::THandle< ND::TSingleHit > > negativeBins, double ChargeExtrapolated);
  virtual ~TTPCHitPad();

  /// Simple analytic fit of the waveform peak.
  /// It uses N bins on each side of the peak
  /// with N = tpcRecon.Reco.Wave.AnalyticFitRange
  void AnalyticFit();

  /// More advanced fitting of the peak using a Gauss
  /// fit from ROOT.
  void RootGaussFit();

  /// Number of peaks in the waveform that could not be split
  /// into individuals waveforms because they overlap too much.
  unsigned int GetNumberPeaks();

  /// Returns a vector with the charge of each one the peaks in the waveform.
  /// If there is only one peak, the given charge is the same as GetCharge().
  std::vector<double> GetPeakCharges();

  /// Returns a vector with the time of each one the peaks in the waveform.
  /// If there is only one peak, the given time is the same as GetTime().
  std::vector<double> GetPeakTimes();

  /// Return a vector with the charges of each of the negative peaks in the waveform
  std::vector<double> GetNegativePeakCharges();
  /// Return a vector with the times of each of the negative peaks in the waveform
  std::vector<double> GetNegativePeakTimes();

// TODO ?
//  void GetWaveformTimeRange(double &min, double &max);

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

  /// Get whether hit is candidate for a hair
  bool GetHairCandidate(){ return fHairCandidate; }

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
  double Y(){return this->GetPosition().Y();};
  /// Short cut to access the Z position of the pad
  double Z(){return this->GetPosition().Z();};

  /// Set whether hit is candidate for a hair
  void SetHairCandidate(bool hairCandidate){ fHairCandidate = hairCandidate; }

  /// Short cut to get the TMultiHit directly
  ND::THandle<ND::THit> GetMultiHit();

  ND::THandle<ND::TMultiHit> ConvertToOAEvent();

private:
  /// Initialise member variables
  void Init();
  void InitParameters();

  /// Vector with all the peak times when the waveform is composed of multipeaks too close to be split up
  std::vector< ND::THandle< ND::TSingleHit > > fPeakBin;
  /// Vector with all the negative peaks in the waveform
  std::vector< ND::THandle< ND::TSingleHit > > fNegativePeakBin;
  /// Index of the bin of maximum charge in the waveform
  int fMaxPeak;
  unsigned int fSaturation;
  int fDriftSense;

  bool fAtHoriEdge;
  bool fAtVertEdge;

  /// Whether hit is candidate for a hair
  bool fHairCandidate;

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

  ClassDef(TTPCHitPad,1);

  //protected:

};

#endif
