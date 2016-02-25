#include "TTPCWaveformDigest.hxx"
#include <TComboHit.hxx>

//*****************************************************************************************************************
TTPCWaveformDigest::TTPCWaveformDigest(void) {
//*****************************************************************************************************************
  fSatThreshold = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.Wave.SaturationThreshold");

  fThresholdSwitch = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.Wave.ThresholdSwitchLevel");
  fRiseThresLowAmp = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.Wave.RiseThresholdLowAmp");
  fFallThresLowAmp = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.Wave.FallThresholdLowAmp");
  fRiseThresHighAmp = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.Wave.RiseThresholdHighAmp");
  fFallThresHighAmp = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.Wave.FallThresholdHighAmp");

  fPeakAmpCut = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.Wave.LDSAmplitudeCut");

  fAnalyticFitRange = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.Wave.AnalyticFitRange");

  fSamplingTime = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.AfterTPC.SamplingTime");

  fPeakSepThreshold = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.Wave.PeakSepThreshold");

  fExtrapolationParameterA = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.Wave.ExtrapolationParameterA");
  fExtrapolationParameterB1 = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.Wave.ExtrapolationParameterB1");
  fExtrapolationParameterB2 = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.Wave.ExtrapolationParameterB2");
}


//*****************************************************************************************************************
TTPCWaveformDigest::~TTPCWaveformDigest() {
//*****************************************************************************************************************
}



//*****************************************************************************************************************
void TTPCWaveformDigest::Process(ND::TND280Event& event, ND::THandle<ND::THitSelection> rhit) {
//*****************************************************************************************************************

  ND::THandle<ND::THitSelection> tpc(event.GetHitSelection("tpc"));
  Process(tpc, rhit);
}



//*****************************************************************************************************************
void TTPCWaveformDigest::Process(const ND::TAlgorithmResult& in, ND::THandle<ND::THitSelection> rhit) {
//*****************************************************************************************************************

  ND::THandle<ND::THitSelection> tpc(in.GetHitSelection("tpc"));
  Process(tpc, rhit);
}


//*****************************************************************************************************************
void TTPCWaveformDigest::Process(ND::THandle<ND::THitSelection> rawtpc, ND::THandle<ND::THitSelection> rhittpc) {
//*****************************************************************************************************************

  if (!rawtpc)
    return;

  for (ND::THitSelection::const_iterator mhit = rawtpc->begin(); mhit != rawtpc->end(); ++mhit) {

    ND::THandle<ND::TMultiHit> mh1 = *mhit;
    ND::THandle<ND::TComboHit> ch1 = *mhit;

    if (mh1){
      //ND::TTPCWaveForm tpcWaveform = LocalDerivativeSearch(mh1);
      // This waveform is way too small, ignore it.
      if (mh1->size() < unsigned(fAnalyticFitRange+1))
        continue;
      std::vector<TTRExSinglePeak> PeakList = MultiMaximumSearch(mh1);
      std::vector<ND::THandle< ND::THit > > tpcWaveforms = BreakWaveforms(mh1, PeakList);

      for (std::vector<ND::THandle< ND::THit > >::iterator wf = tpcWaveforms.begin(); wf != tpcWaveforms.end(); wf++) {
        rhittpc->push_back(*wf);
      }
    }
    else if (ch1){
      ND::THandle<ND::TComboHit> newcombo(new ND::TComboHit);
      for (ND::THitSelection::const_iterator hpad = ch1->GetHits().begin(); hpad != ch1->GetHits().end();hpad++){
        ND::THandle<ND::TMultiHit> mh2 = *hpad;
        // This waveform is way too small, ignore it.
        if (mh2){
          if (mh2->size() < unsigned(fAnalyticFitRange+1))
            continue;
          //ND::TTPCWaveForm tpcWaveform = LocalDerivativeSearch(mh2);
          std::vector<TTRExSinglePeak> PeakList = MultiMaximumSearch(mh2);
          std::vector<ND::THandle< ND::THit > > tpcWaveforms = BreakWaveforms(mh2, PeakList);

          for (std::vector<ND::THandle< ND::THit > >::iterator wf = tpcWaveforms.begin(); wf != tpcWaveforms.end(); wf++) {
            newcombo->AddHit(*wf);
          }
        }
        else {
          ND280Severe("TTPCWaveformDigest: input TPC hit is neither a TMultiHit, nor a TComboHit.");
          throw ND::ETRExHitNotTMultiHit();
        }
      }
      rhittpc->push_back(newcombo);
    }
    else {
      ND280Severe("TTPCWaveformDigest: input TPC hit is neither a TMultiHit, nor a TComboHit.");
      throw ND::ETRExHitNotTMultiHit();
    }

  }

}


//*****************************************************************************************************************
void TTPCWaveformDigest::Process(ND::THandle<ND::THitSelection> rawtpc, ND::THitSelection *rhittpc) {
//*****************************************************************************************************************

  if (!rawtpc)
    return;

  for (ND::THitSelection::const_iterator mhit = rawtpc->begin(); mhit != rawtpc->end(); ++mhit) {

    ND::THandle<ND::TMultiHit> mh1 = *mhit;
    ND::THandle<ND::TComboHit> ch1 = *mhit;

    if (mh1){
      //ND::TTPCWaveForm tpcWaveform = LocalDerivativeSearch(mh1);
      // This waveform is way too small, ignore it.
      if (mh1->size() < unsigned(fAnalyticFitRange+1))
        continue;
      std::vector<TTRExSinglePeak> PeakList = MultiMaximumSearch(mh1);
      std::vector<ND::THandle< ND::THit > > tpcWaveforms = BreakWaveforms(mh1, PeakList);

      for (std::vector<ND::THandle< ND::THit > >::iterator wf = tpcWaveforms.begin(); wf != tpcWaveforms.end(); wf++) {
        rhittpc->push_back(*wf);
      }
    }
    else if (ch1){
      ND::THandle<ND::TComboHit> newcombo(new ND::TComboHit);
      for (ND::THitSelection::const_iterator hpad = ch1->GetHits().begin(); hpad != ch1->GetHits().end();hpad++){
        ND::THandle<ND::TMultiHit> mh2 = *hpad;
        // This waveform is way too small, ignore it.
        if (mh2){
          if (mh2->size() < unsigned(fAnalyticFitRange+1))
            continue;
          //ND::TTPCWaveForm tpcWaveform = LocalDerivativeSearch(mh2);
          std::vector<TTRExSinglePeak> PeakList = MultiMaximumSearch(mh2);
          std::vector<ND::THandle< ND::THit > > tpcWaveforms = BreakWaveforms(mh2, PeakList);

          for (std::vector<ND::THandle< ND::THit > >::iterator wf = tpcWaveforms.begin(); wf != tpcWaveforms.end(); wf++) {
            newcombo->AddHit(*wf);
          }
        }
        else {
          ND280Severe("TTPCWaveformDigest: input TPC hit is neither a TMultiHit, nor a TComboHit.");
          throw ND::ETRExHitNotTMultiHit();
        }
      }
      rhittpc->push_back(newcombo);
    }
    else {
      ND280Severe("TTPCWaveformDigest: input TPC hit is neither a TMultiHit, nor a TComboHit.");
      throw ND::ETRExHitNotTMultiHit();
    }

  }

}

//*****************************************************************************************************************
ND::THandle<ND::THitSelection> TTPCWaveformDigest::Process(ND::THandle<ND::THitSelection> rawtpc) {
//*****************************************************************************************************************
  ND::THandle<ND::THitSelection> rhits(new ND::THitSelection());
  Process(rawtpc, rhits);
  return rhits;
}

//*****************************************************************************************************************
std::vector<TTRExSinglePeak> TTPCWaveformDigest::LocalDerivativeSearch(ND::THandle<ND::TMultiHit> mh) {
//*****************************************************************************************************************

  int numPeaks = 0;

  std::vector<TTRExSinglePeak> PeaksFound;

  typedef std::vector< ND::THandle< ND::TSingleHit > >::const_iterator iterator;

  double charge, time;
  int i = 1;
  // Skip the end time bins as alogrithm requires there to always be defined
  // one lower and higher time bin then any bin being examined for as a potential
  // peak
  for (iterator bin = mh->begin()+1; bin !=  mh->end()-1; bin++) {

    charge = (*bin)->GetCharge();
    time = (*bin)->GetTime();

    // Cut on the bin amplitude to avoid finding non-existent peaks
    // in the pedestal time bins
    if (charge < fPeakAmpCut) {
      i++;
      continue;
    }

    // Look at the local derivative of the waveform. If the current
    // bin is larger than the four nearest bins then a peak has been
    // found and is registered in lists:
    if ((*(bin-1))->GetCharge() < charge && (*(bin+1))->GetCharge() < charge) {
      if ((*(bin-2))->GetCharge() < charge && (*(bin+2))->GetCharge() < charge) {
        TTRExSinglePeak Peak(charge, i, time, 0, charge);
        PeaksFound.push_back(Peak);
        numPeaks++;
        bin += 2;
        i += 2;
      }
    }
    // Failing the above check look for a saturated peak  where
    // multiple bins have equal amplitude.
    else if (charge == (*(bin+1))->GetCharge()) {

      if ((*(bin-1))->GetCharge() < charge && (*(bin-2))->GetCharge() < charge) {

        int j = i;
        iterator bin2 = bin;
        do {
          if (bin2 == mh->end()) {
            break;
          }
          bin2 += 1;
          j++;
        } while (charge == (*bin2)->GetCharge());

        if ((*bin2)->GetCharge() < charge) {

          if (charge < fSatThreshold) {
            TTRExSinglePeak Peak(charge, i, time, 0, charge);
            PeaksFound.push_back(Peak);
          }
          else {
            int index = int(ceil((j-i)/2.0)) - 1;
            int t = int(time + fSamplingTime*(ceil((j-i)/2.0) - 1));
            index += i;

            // Charge extrapolation for saturated peaks
            double CF = 1.;
            if (bin2 != (mh->end()-1)){ // Is there a bin right of bin2?
              CF = SaturationCorrectionFactor(j-i, (*(bin-2))->GetCharge(), (*(bin-1))->GetCharge(), (*(bin2))->GetCharge(), (*(bin2+1))->GetCharge() );
            } else { // No there is not.
              CF = SaturationCorrectionFactor(j-i, (*(bin-2))->GetCharge(), (*(bin-1))->GetCharge(), (*(bin2))->GetCharge(), 0. );
            }

            // The extrapolation is based on the higher of the two neighbouring bins.
            double echarge;
            if ((*(bin-1))->GetCharge() > (*(bin2))->GetCharge()){
              echarge = (*(bin-1))->GetCharge() * CF;
            } else {
              echarge = (*(bin2))->GetCharge() * CF;
            }
            // Use extrapolated charge only if it is higher than the saturated charge.
            if (echarge < charge){
              echarge = charge;
            }

            TTRExSinglePeak Peak(charge, index, t, j-i, echarge);
            PeaksFound.push_back(Peak);
          }

          numPeaks++;
          i = j;
          bin = bin2;
        }
      }
    }

    i++;
  }

  // If the algorithm failed to find any peaks, run the old single peak
  // search code that looks for the first bin with the largest charge
  if (numPeaks == 0) {
    std::cout << "Calling fallback..." << std::endl;

    TTRExSinglePeak Peak(-99999.9, -1, -1, 0, -99999.9);
    if (FindStrongestPeak(mh, Peak)) {
      PeaksFound.push_back(Peak);
      numPeaks++;
    }
  }

  return PeaksFound;
}

//*****************************************************************************************************************
std::vector<TTRExSinglePeak> TTPCWaveformDigest::MultiMaximumSearch(ND::THandle<ND::TMultiHit> mh, bool findNegative) {
//*****************************************************************************************************************
  int numPeaks = 0;

  std::vector<TTRExSinglePeak> PeaksFound;

  typedef std::vector< ND::THandle< ND::TSingleHit > >::const_iterator iterator;

  // Variables used in determining when one peak is found an new search
  // should be started they are the number of consecutive waveform bins
  // where the slope was increasing or decreasing
  int increasing = 0,
    decreasing = 0;

  // Add comment on variables.
  double maxCharge = -99999.9,
    prevCharge = 0,
    maxTime = 0;
  int maxBin = 0,
    satCount = 0,
    index = 1;
  bool findPeak = false;

  if(findNegative){
    maxCharge = 99999.9;
  };

  prevCharge = (*(mh->begin()))->GetCharge();
  for (iterator bin = mh->begin()+1; bin != mh->end(); bin++) {

    double charge = (*bin)->GetCharge();
    double time = (*bin)->GetTime();

    if (findPeak) {
      if ( (!findNegative && (charge > maxCharge)) || (findNegative && (charge < maxCharge)) ) {
        maxCharge = charge;
        maxTime = time;
        maxBin = index;
        decreasing = 0;
        satCount = 0;
      }
      else if (charge == maxCharge) {
        if (satCount == 0) {
          satCount = 1;
        }
        satCount++;
        decreasing = 0;
      }
      else {
        decreasing++;

        // Wavefrom has decreased an amount such that the maximum can be
        // taken as a found peak!
        bool decreasedEnough;
        if(findNegative){
          // no idea how to work out what the switch should be
          decreasedEnough = (decreasing >= fFallThresHighAmp);
        }
        else{
          decreasedEnough = ((charge < fThresholdSwitch && decreasing >= fFallThresLowAmp) ||
            (charge > fThresholdSwitch && decreasing >= fFallThresHighAmp));
        };
        if (decreasedEnough) {
          double echarge = maxCharge;
          if (!findNegative && maxCharge < fSatThreshold) {
            satCount = 0;
          }
          else if (satCount != 0) {
            if(!findNegative){
              // Charge extrapolation for saturated peaks
              double l2, l1, r1, r2;
              if(maxBin >= 2){
                l2 = (*(mh->begin()+maxBin-2))->GetCharge();
              }else{
                l2 = 0;
              }
              if(maxBin >= 1){
                l1 = (*(mh->begin()+maxBin-1))->GetCharge();
              }else{
                l1 = 0;
              }
              if(maxBin+satCount <= int(mh->size()-1)){
                r1 = (*(mh->begin()+maxBin+satCount))->GetCharge();
              }else{
                r1 = 0;
              }
              if(maxBin+satCount+1 <= int(mh->size()-1)){
                r2 = (*(mh->begin()+maxBin+satCount+1))->GetCharge();
              }else{
                r2 = 0;
              }
              double CF = SaturationCorrectionFactor(satCount, l2, l1, r1, r2);
              if (l1 > r1){
                echarge = l1 * CF;
              } else {
                echarge = r1 * CF;
              }
              if (echarge < maxCharge){
                echarge = maxCharge;
              }
            };

            int shift = int(ceil(satCount/2.0)) - 1;
            maxTime += fSamplingTime*shift;
            maxBin += shift;
          };
          if( (!findNegative && (maxCharge > 0.)) || (findNegative && (maxCharge < 0.)) ){
            TTRExSinglePeak Peak(maxCharge, maxBin, maxTime, satCount, echarge);
            PeaksFound.push_back(Peak);
            numPeaks++;
          };

          findPeak = false;
          if(findNegative){
            maxCharge = 99999.9;
          }
          else{
            maxCharge = -99999.9;
          };

          decreasing = 0;
          satCount = 0;
        };
      };
    }
    // Look for the start of a new peak, indicated by an upward turn in the slope.
    else {
      if ( (!findNegative && (charge > prevCharge)) || (findNegative && (charge < prevCharge)) ) {
        increasing++;

        // Once the threshold is reached for a steady increase in the signal
        // start a new "peak search"
        bool increasedEnough;
        if(findNegative){
          // no idea how to work out what the switch should be
          increasedEnough = (increasing >= fRiseThresHighAmp);
        }
        else{
          increasedEnough = ((charge < fThresholdSwitch && increasing >= fRiseThresLowAmp) ||
            (charge > fThresholdSwitch && increasing >= fRiseThresHighAmp));
        };
        if (increasedEnough) {
          increasing = 0;
          prevCharge = 0;
          findPeak = true;

          // This step prevents the algorithm from missing a peak because once it decides
          // to start searching for another peak it is already at the peak
          bin--;
          index--;
        }
      }
      else {
        // Requires the increase in amplitude be continuous before new "peak search"
        increasing = 0;
      }
    }

    index++;
    prevCharge = charge;
  }

  // If the algorithm failed to find any peaks, run the old single peak
  // search code that looks for the first bin with the largest charge
  if (numPeaks == 0) {

    TTRExSinglePeak Peak(-99999.9, -1, -1, 0, -99999.9);
    if (FindStrongestPeak(mh, Peak, findNegative)) {
      PeaksFound.push_back(Peak);
      numPeaks++;
    }
  }

  return PeaksFound;
}


//*****************************************************************************************************************
bool TTPCWaveformDigest::FindStrongestPeak(ND::THandle<ND::TMultiHit> mh, TTRExSinglePeak &peak, bool findNegative) {
//*****************************************************************************************************************

  typedef std::vector< ND::THandle< ND::TSingleHit > > Hits;
  typedef Hits::const_iterator iterator;

  double maxCharge = -99999.9;
  double l2Charge = 0.0;
  double l1Charge = 0.0;
  double r1Charge = 0.0;
  double r2Charge = 0.0;
  int maxTime = 0,
    maxBin = 0,
    satCount = 0,
    index = 0;

  if(findNegative){
    maxCharge = 99999.9;
  };

  for (iterator bin = mh->begin(); bin != mh->end(); bin++) {
    double charge = (*bin)->GetCharge();
    double time = (*bin)->GetTime();

    if ( (!findNegative && (charge > maxCharge)) || (findNegative && (charge < maxCharge)) ) {
      l2Charge = l1Charge;
      if ( (!findNegative && (maxCharge < 0.)) || (findNegative && (maxCharge > 0.)) ){ // Prevent default value from disturbing the algorithm
        l1Charge = 0.;
      } else {
        l1Charge = maxCharge;
      }
      maxCharge = charge;
      maxTime = int(time);
      maxBin = index;
      satCount = 0;
      r1Charge = 0.0;
      r2Charge = 0.0;
    }
    // Code for peak saturation:
    else if (charge == maxCharge) {
      if (satCount == 0) {
        satCount = 1;
      }
      satCount++;
    }
    // Also get the two bins right of the peak
    else {
      if (satCount > 0){
        if (index == maxBin + satCount){
          r1Charge = charge;
        } else if (index == maxBin + satCount + 1){
          r2Charge = charge;
        }
      } else {
        if (index == maxBin + 1){
          r1Charge = charge;
        } else if (index == maxBin + 2){
          r2Charge = charge;
        }
      }
    }

    index++;
  }

  // Store the peak location:
  if ( (!findNegative && (maxCharge > -9999.9)) || (findNegative && (maxCharge < 9999.9)) ) {
    double peakCharge,
           peakTime,
           peakExtrapolatedCharge;
    int peakBin,
        peakSatCount;
    peakCharge = maxCharge;

    // If not above a threshold then peak is not saturated but just has
    // multiple bins of equal amplitude
    if (!findNegative && peakCharge < fSatThreshold) {
      peakSatCount = 0;
      peakTime = maxTime;
      peakBin = maxBin;
      peakExtrapolatedCharge = peakCharge;
    }
    else {
      int shift = int(ceil(satCount/2.0)) - 1;
      peakSatCount = satCount;
      peakTime = int(maxTime + fSamplingTime*shift);
      peakBin = maxBin + shift;

      if(findNegative){
        peakExtrapolatedCharge = peakCharge;
      }
      else{
        // Charge extrapolation for saturated peaks
        double CF = SaturationCorrectionFactor( peakSatCount, l2Charge, l1Charge, r1Charge, r2Charge );
        // The extrapolation is based on the higher of the two neighbouring bins.
        double echarge;
        if (l1Charge > r1Charge){
          echarge = l1Charge * CF;
        } else {
          echarge = r1Charge * CF;
        }
        // Use extrapolated charge if it is higher than the saturated charge.
        if (echarge > peakCharge){
          peakExtrapolatedCharge = echarge;
        } else {
          peakExtrapolatedCharge = peakCharge;
        }
      }
    }

    // Something went wrong. Probably a messed up waveform found in data.
    // Just make sure that things don't crash and this waveform won't be used downstream anyway.
    if (peakBin < 0){
      peakBin = 0;
    } else if (peakBin >= index){
      peakBin = index -1;
    }
    peak.InitPeak(peakCharge, peakBin, peakTime, peakSatCount, peakExtrapolatedCharge);
  }
  // This should never happen, unless the waveform has less than fAnalyticFitRange bins...
  else {
    std::cerr << "TTPCWaveformDigest::FindStrongestPeak(): Warning!! Fallback peak search found 0 peaks (not possible?)!" << std::endl;
    return false;
  }

  return true;
}

//*****************************************************************************************************************
int TTPCWaveformDigest::FindMinimum(ND::THandle<ND::TMultiHit> mh, TTRExSinglePeak low, TTRExSinglePeak high, bool findNegative) {
//*****************************************************************************************************************

  typedef std::vector< ND::THandle< ND::TSingleHit > > Hits;
  typedef Hits::const_iterator iterator;

  double minCharge = 99999.9;
  int minBin = 0,
    index,
    start = low.Bin(),
    end = high.Bin();
  index = start;

  if(findNegative){
    minCharge = -99999.9;
  };

  for (iterator bin = mh->begin() + start; bin != mh->begin() + end; bin++) {
    double charge = (*bin)->GetCharge();

    if ( (!findNegative && (charge < minCharge)) || (findNegative && (charge > minCharge)) ) {
      minCharge = charge;
      minBin = index;
    }
    index++;
  }

  // Checks that there is a true minimum, such that the time bin with the lowest
  // amplitude is not one of those at the ends of the range bounded by the provided
  // peaks.
  if (minBin == start) {
    return -1;
  }
  else if (minBin == end) {
    return -2;
  }

  return minBin;
}

//*****************************************************************************************************************
std::vector< ND::THandle< ND::THit > > TTPCWaveformDigest::BreakWaveforms( ND::THandle<ND::TMultiHit> mh, std::vector<TTRExSinglePeak> Peaks) {
//*****************************************************************************************************************
  typedef std::vector< ND::THandle< ND::TSingleHit > > Bins;
  typedef Bins::const_iterator iterator;

  std::vector< ND::THandle< ND::THit > > waveforms;
  if (Peaks.size() < 1)
    return waveforms;

  Bins closePeaks;
  double closePeakMaximum = -0xABCDEF;
  double closePeakMaximumExtrapolated = -0xABCDEF;
  double closePeakMaximumSaturation = 0;

  // One peak, straight forward
  if (Peaks.size() == 1) {
    ND::THandle<ND::THit> hit = mh;
    ND::TWritableReconHit whit(hit);

    Bins minSinglePeaks = GetNegatives(mh, Peaks[0]);

    whit.SetCharge(Peaks[0].Charge());
    whit.SetTime(Peaks[0].Time());
    whit.SetTimeUncertainty( Peaks[0].Saturation() * fSamplingTime /2.);
    TVector3 tmpPos = hit->GetPosition();
    TVector3 tmpUnc = hit->GetUncertainty();
    whit.SetPosition(tmpPos);
    whit.SetUncertainty(tmpUnc);
    closePeaks.push_back((*mh)[Peaks[0].Bin()]);
    ND::THandle<ND::TTPCHitPad> rhit(new ND::TTPCHitPad(whit, closePeaks, minSinglePeaks, Peaks[0].ExtrapolatedCharge()));
    //ND::THandle<ND::TTPCHitPad> rhit(new ND::TTPCHitPad(whit, closePeaks, Peaks[0].ExtrapolatedCharge()));
    waveforms.push_back(rhit);
  }
  // Make a new TMultiHit for each peak with corresponding bins copied
  else {
    iterator firstBin = mh->begin();
    iterator lastBin = mh->begin();
    for (unsigned int p = 0; p < (Peaks.size()-1); p++) {
      TTRExSinglePeak thisP = Peaks[p]; TTRExSinglePeak nextP = Peaks[p+1];

      // The peaks are too close, don't split them
      if (( (nextP.Bin() - thisP.Bin()) + nextP.Saturation() + thisP.Saturation() ) < fPeakSepThreshold){
        // Extract the TSingleHit corresponding to this peak
        closePeaks.push_back((*mh)[thisP.Bin()]);
        // Save the highest peak charge
        if (thisP.ExtrapolatedCharge() > closePeakMaximumExtrapolated){
          closePeakMaximum = thisP.Charge();
          closePeakMaximumExtrapolated = thisP.ExtrapolatedCharge();
          closePeakMaximumSaturation = thisP.Saturation();
        }
        continue;
      }

      int minBin = FindMinimum(mh, thisP, nextP);

      if (minBin > 0){
        lastBin = mh->begin() + minBin + 1;
      }
      else {      // Use the mid point between the two peaks as cut off since we didn't find a minimum
        double timeCut = thisP.Time() + ((nextP.Time() - thisP.Time()) / 2.);
        // stop the waveform at the last bin _before_ the cut off.
        for (; lastBin != mh->end(); lastBin++) {
          if ( (*lastBin)->GetTime() > timeCut) {
             break;
          }
        }
      }
      ND::THandle<ND::TMultiHit> newHit(new ND::TMultiHit(firstBin, lastBin));

      // The TTPCHitPad constituent is the broken waveform
      ND::THandle<ND::THit> hit = newHit;
      ND::TWritableReconHit whit(hit);
      whit.AddHit(mh);

      // Save the TSingleHits corresponding to each one of the close peaks in this waveform
      double waveformTime = Peaks[p].Time();
      double waveformSaturation = Peaks[p].Saturation();
      double waveformCharge = Peaks[p].Charge();
      double waveformChargeExtrapolated = Peaks[p].ExtrapolatedCharge();

      if (closePeaks.size() > 0){
        // The last peak hasn't been added to the list yet
        closePeaks.push_back((*mh)[thisP.Bin()]);
        // Save the highest peak charge
        if (thisP.ExtrapolatedCharge() > closePeakMaximumExtrapolated){
          closePeakMaximum = thisP.Charge();
          closePeakMaximumExtrapolated = thisP.ExtrapolatedCharge();
          closePeakMaximumSaturation = thisP.Saturation();
        }
        double firstTime = closePeaks[0]->GetTime();
        double lastTime = closePeaks[closePeaks.size()-1]->GetTime();
        // Many close peaks so take the mid point between the first and last peaks
        waveformTime = (lastTime + firstTime)/2.;
        // for the rest use the values from the highest peak
        waveformSaturation = closePeakMaximumSaturation;
        waveformCharge = closePeakMaximum;
        waveformChargeExtrapolated = closePeakMaximumExtrapolated;
      }

      Bins minSinglePeaks = GetNegatives(newHit, thisP);

      whit.SetCharge(waveformCharge);
      whit.SetTime(waveformTime);
      whit.SetTimeUncertainty( waveformSaturation * fSamplingTime /2.);
      TVector3 tmpPos = hit->GetPosition();
      TVector3 tmpUnc = hit->GetUncertainty();
      whit.SetPosition(tmpPos);
      whit.SetUncertainty(tmpUnc);
      ND::THandle<ND::TTPCHitPad> rhit(new ND::TTPCHitPad(whit, closePeaks, minSinglePeaks, waveformChargeExtrapolated));
      waveforms.push_back(rhit);

      // Reset close peak finder
      closePeakMaximum = -0xABCDEF;
      closePeakMaximumExtrapolated = -0xABCDEF;
      closePeakMaximumSaturation = 0;
      closePeaks.clear();
      firstBin = lastBin;
    }
    // Take care of the last peak
    lastBin = mh->end();
    int p = Peaks.size()-1;

    ND::THandle<ND::TMultiHit> newHit(new ND::TMultiHit(firstBin, lastBin));
    // The TTPCHitPad constituent is the broken waveform
    ND::THandle<ND::THit> hit = newHit;
    ND::TWritableReconHit whit(hit);
    whit.AddHit(mh);

    // Save the TSingleHits corresponding to each one of the close peaks in this waveform
    double waveformTime = Peaks[p].Time();
    double waveformSaturation = Peaks[p].Saturation();
    double waveformCharge = Peaks[p].Charge();
    double waveformChargeExtrapolated = Peaks[p].ExtrapolatedCharge();

    if (closePeaks.size() > 0){
      // The last peak hasn't been added to the list yet
      closePeaks.push_back((*mh)[Peaks[p].Bin()]);
      // Save the highest peak charge
      if (Peaks[p].ExtrapolatedCharge() > closePeakMaximumExtrapolated){
        closePeakMaximum = Peaks[p].Charge();
        closePeakMaximumExtrapolated = Peaks[p].ExtrapolatedCharge();
        closePeakMaximumSaturation = Peaks[p].Saturation();
      }
      double firstTime = closePeaks[0]->GetTime();
      double lastTime = closePeaks[closePeaks.size()-1]->GetTime();
      // Many close peaks so take the mid point between the first and last peaks
      waveformTime = (lastTime + firstTime)/2.;
      // for the rest use the values from the highest peak
      waveformSaturation = closePeakMaximumSaturation;
      waveformCharge = closePeakMaximum;
      waveformChargeExtrapolated = closePeakMaximumExtrapolated;
    }

    Bins minSinglePeaks = GetNegatives(mh, Peaks[p]);

    whit.SetCharge(waveformCharge);
    whit.SetTime(waveformTime);
    whit.SetTimeUncertainty( waveformSaturation * fSamplingTime /2.);
    TVector3 tmpPos = hit->GetPosition();
    TVector3 tmpUnc = hit->GetUncertainty();
    whit.SetPosition(tmpPos);
    whit.SetUncertainty(tmpUnc);
    ND::THandle<ND::TTPCHitPad> rhit(new ND::TTPCHitPad(whit, closePeaks, minSinglePeaks, waveformChargeExtrapolated));
    waveforms.push_back(rhit);
  }

  return waveforms;
}

//*****************************************************************************************************************
double TTPCWaveformDigest::SaturationCorrectionFactor( int nsat, double left2, double left1, double right1, double right2 ) {
//*****************************************************************************************************************
  // Paramters for the extrapolation
  double parA = fExtrapolationParameterA; // 1.24e-2;
  double parB1 = fExtrapolationParameterB1; // -8.36e-2;
  double parB2 = fExtrapolationParameterB2; // 9.46e-2;

  // The correction factor
  double CF = 1.;

  // The slopes left and right of the plateau
  double left_slope = left1 - left2;
  double right_slope = right2 - right1;

  // Define the plateau as the higher of the two neighbouring bins,
  // since we cannot trust the saturated peaks.
  double plateau_width = nsat + 1;
  double plateau_width_correction = 0;
  if (left1 > right1){ // Left bin defines the plateau.
    if (right_slope < 0.){ // Is the slope reasonable?
      plateau_width_correction = (right1 - left1) / right_slope;
    }
  } else { // Right bin defines the plateau.
    if (left_slope > 0.){ // Is the slope reasonable?
      plateau_width_correction = (right1 - left1) / left_slope;
    }
  }

  // First step correction factor calculation
  CF = exp(parA * pow(plateau_width - plateau_width_correction, 2));

  // Second step correction factor correction
  CF *= 1 + parB1 + parB2*pow(plateau_width_correction, 2);

  // Make sure the correction factor is at all reasonable
  if (CF < 1.){
    CF = 1.;
  }

  return CF;
}

std::vector< ND::THandle<ND::TSingleHit> > TTPCWaveformDigest::GetNegatives(ND::THandle<ND::TMultiHit> mh, TTRExSinglePeak peak){
  std::vector< ND::THandle<ND::TSingleHit> > minPeaks;
  std::vector< ND::THandle<ND::TMultiHit> > checkHits;
  // explicitly search either side of the peak
  if(peak.Bin() > 0 && peak.Bin() < (int)mh->size()-1){
    ND::THandle<ND::TMultiHit> startHit(new ND::TMultiHit(mh->begin(), mh->begin()+peak.Bin()));
    ND::THandle<ND::TMultiHit> endHit(new ND::TMultiHit(mh->begin()+peak.Bin()+1, mh->end()));

    checkHits.push_back(startHit);
    checkHits.push_back(endHit);
  }
  else{
    ND::THandle<ND::TMultiHit> checkHit(new ND::TMultiHit(mh->begin(), mh->end()));

    checkHits.push_back(checkHit);
  };

  for(std::vector< ND::THandle<ND::TMultiHit> >::iterator checkHitIt = checkHits.begin(); checkHitIt != checkHits.end(); ++checkHitIt){
    ND::THandle<ND::TMultiHit> checkHit = *checkHitIt;
    if((int)checkHit->size() > fAnalyticFitRange){
      std::vector<TTRExSinglePeak> minTPCPeaks = MultiMaximumSearch(checkHit, true);
      for(std::vector<TTRExSinglePeak>::iterator minTPCPeakIt = minTPCPeaks.begin(); minTPCPeakIt != minTPCPeaks.end(); ++minTPCPeakIt){
        minPeaks.push_back((*checkHit)[minTPCPeakIt->Bin()]);
      };
    };
  };

  return minPeaks;
}
