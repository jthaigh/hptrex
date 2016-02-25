#ifndef TTPCWaveformDigest_hxx_seen
#define TTPCWaveformDigest_hxx_seen

#include <HEPUnits.hxx>
#include <THandle.hxx>
#include "TRealDatum.hxx"

#include "TOADatabase.hxx"
#include <TComboHit.hxx>
#include <THandle.hxx>
#include <THit.hxx>
#include <TMultiHit.hxx>
#include <TAlgorithm.hxx>
#include <THandle.hxx>
#include <TMultiHit.hxx>
#include <TND280Event.hxx>
#include <TOARuntimeParameters.hxx>
#include <TTPCDigit.hxx>

#include "TTPCHitPad.hxx"

#include <list>
#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <sstream>
#include <stdio.h>


namespace ND {
  // If for some reason an input TPC hit is not a TMultiHit
  OA_EXCEPTION(ETRExHitNotTMultiHit,EoaCore);
}
;

class TTRExSinglePeak;

class TTPCWaveformDigest {
public:
  TTPCWaveformDigest();

  virtual ~TTPCWaveformDigest();

  void Process(const ND::TAlgorithmResult& in, ND::THandle<ND::THitSelection> rhit);

  void Process(ND::TND280Event& event, ND::THandle<ND::THitSelection> rhit);

  void Process( ND::THandle<ND::THitSelection> rawtpc, ND::THandle<ND::THitSelection> rhittpc );
  void Process( ND::THandle<ND::THitSelection> rawtpc, ND::THitSelection *rhittpc );

  ND::THandle<ND::THitSelection> Process(ND::THandle<ND::THitSelection> rawtpc);


  /**
   * Performs an iterative search over the time bins of the waveform stored in the
   * TMultiHit looking for the signal peaks.
   *
   * Method breaks the wavefrom up into sections, based on finding the rising and
   * falling edges of signal peaks. Within reach section, or peak region, the algorithm
   * then looks to find the time bin with the largest amplitude. In the event that
   * multiple time bins are of equal amplitude the algorithm tracks the number of
   * saturated bins if:
   *
   *                           charge > fSatThreshold
   *
   * the peak is flagged as saturated and peak position is revised, from the first
   * of the saturated bins, to the median bin of the peak region, if there is an
   * odd number of bins, or the floored average for an even number of bins.
   *
   * In the event that no peak is found then the single peak search is called.
   *
   */
  std::vector<TTRExSinglePeak > MultiMaximumSearch( ND::THandle<ND::TMultiHit> mh, bool findNegative=false);


  /*, bool negativeCharge=false*
   * The original waveform is "broken" into waveforms containing the bins
   * corresponding to each peak found in the peak search.
   *
   */
  std::vector<ND::THandle< ND::THit > > BreakWaveforms( ND::THandle<ND::TMultiHit> mh, std::vector<TTRExSinglePeak> Peaks);


  /**
   * Calculates a correction factor for the charge of saturated peaks,
   * using the number of saturated bins and two bins each left and right of the plateau.
   *
   */
  double SaturationCorrectionFactor( int nsat, double left2, double left1, double right1, double right2 );


  /**
   * Performs an iterative search over the time bins of the waveform stored in the
   * TMultiHit looking for the signal peaks.
   *
   * This method looks at the waveform time bins with amplitude:
   *
   *                           charge > fPeakAmpCut
   *
   * and uses the local derivative to identify peaks. The bin's charge is compared to
   * those of the two neighbouring bins on either side. If the time bin has a larger
   * amplitude then the neighbours then a peak registered. In the event that multiple
   * time bins are of equal amplitude the algorithm tracks the number of saturated
   * bins if:
   *
   *                           charge > fSatThreshold
   *
   * the peak is flagged as saturated and peak position is revised, from the first
   * of the saturated bins, to the median bin of the peak region, if there is an
   * odd number of bins, or the floored average for an even number of bins.
   *
   * In the event that no peak is found then the single peak search is called.
   *
   */
  std::vector<TTRExSinglePeak> LocalDerivativeSearch(ND::THandle<ND::TMultiHit> mh);


  /**
   * Performs an iterative search over the time bins of the waveform stored in the
   * TMultiHit and looks for a single signal peak.
   *
   * The highest charge time bin is found and taken to be the signal peak, assuming that the
   * waveform contains only one peak. In the event that multiple time bins are of equal
   * amplitude the algorithm tracks the number of saturated bins if:
   *
   *                           charge > fSatThreshold
   *
   * The peak then flagged as saturated and peak position is revised, from the first of the
   * saturated bins, to the median bin of the peak region if there is an odd number of
   * bins, or the floored average for an even number of bins.
   *
   * @param mh the TMultiHit the contains the TPC wavefrom
   * @param peak the TTRExSinglePeak that will be filled with the found peak
   */
  bool FindStrongestPeak(ND::THandle<ND::TMultiHit> mh, TTRExSinglePeak &Peak, bool findNegative=false);

  /**
    * Searches for the minimum in the waveform stored the provided TMultiHit.
    *
    * Searches the waveform between the two provided waveform peaks for the time bin with
    * the minimum signal amplitude. As a check that there is minimum between the provided
    * peaks the bin, found to have the lowest signal, is confirmed to not be either of the
    * provided peaks locations. These conditions should not be satisfied for peak
    * positions returned from any of the multiple peak finding functions, but are checked
    * for in the event of errors.
    *
    * Note, this function does not check the separation distance between the minimum and
    * the provided peaks. If used for breaking of the waveform it is left to the code
    * employing this function to check that the minimum is far enought from the given peaks
    * to not interfere with subsquent operations, such as fiting of the peaks .
    *
    * @param mh the TMultiHit the contains the TPC wavefrom
    * @param low peak located at the eariler time
    * @param high peak located at the later time
    * @return time bin index in the waveform for the location of the minimum signal, or -1 if min = low, or -2 if min = high
    * @see #MultiMaximumSearch( ND::THandle<ND::TMultiHit> mh)
    * @see #LocalDerivativeSearch(ND::THandle<ND::TMultiHit> mh)
    */
  int FindMinimum(ND::THandle<ND::TMultiHit> mh, TTRExSinglePeak low, TTRExSinglePeak high, bool findNegative=false);

private:
  /// trexRecon.Reco.Wave.RiseThresholdLowAmp
  int fRiseThresLowAmp,
  /// trexRecon.Reco.Wave.FallThresholdLowAmp
  fFallThresLowAmp,
  /// trexRecon.Reco.Wave.RiseThresholdHighAmp
  fRiseThresHighAmp,
  /// trexRecon.Reco.Wave.FallThresholdHighAmp
  fFallThresHighAmp,
  /// trexRecon.Reco.Wave.AnalyticFitRange
  fAnalyticFitRange;

  /// trexRecon.Reco.Wave.SaturationThreshold
  double fSatThreshold,
  /// trexRecon.Reco.Wave.ThresholdSwitchLevel
  fThresholdSwitch,
  /// trexRecon.Reco.Wave.LDSAmplitudeCut
  fPeakAmpCut;

  /// Sampling time of the TPC waveforms
  double fSamplingTime;

  /// trexRecon.Reco.Wave.PeakSepThreshold
  int fPeakSepThreshold;

  /// trexRecon.Reco.Wave.ExtrapolationParameterA
  double fExtrapolationParameterA;
  /// trexRecon.Reco.Wave.ExtrapolationParameterB1
  double fExtrapolationParameterB1;
  /// trexRecon.Reco.Wave.ExtrapolationParameterB2
  double fExtrapolationParameterB2;

  /// Find negative peaks from a waveform
  std::vector< ND::THandle<ND::TSingleHit> > GetNegatives(ND::THandle<ND::TMultiHit> mh, TTRExSinglePeak peak);

};

/// Selected peak of the waveform.
class TTRExSinglePeak {
  public:
    TTRExSinglePeak(double charge, int bin, double time, int saturation, double extrapolatedCharge){
      InitPeak(charge, bin, time, saturation, extrapolatedCharge);
    };
    void InitPeak(double charge, int bin, double time, int saturation, double extrapolatedCharge){
      fBin = bin;
      fCharge = charge;
      fTime = time;
      fSaturation = saturation;
      fExtrapolatedCharge = extrapolatedCharge;
    };
    int Bin() {return fBin;};
    double Charge() {return fCharge;};
    double Time() {return fTime;};
    int Saturation() {return fSaturation;};
    double ExtrapolatedCharge() {return fExtrapolatedCharge;};

  private:
    int fBin;
    double fCharge;
    double fTime;
    int fSaturation;
    double fExtrapolatedCharge;
};

#endif
