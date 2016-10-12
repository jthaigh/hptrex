#ifndef TRexReco_hxx_seen
#define TRexReco_hxx_seen

#include "patternRecognition/TTPCTRExPatAlgorithm.hxx"

#include "TTPCSeeding.hxx"
#include "TTPCT0Finder.hxx"
#include "TTPCTracking.hxx"
#include "TTPCLikelihoodMatch.hxx"
#include "TTPCPID.hxx"
#include "TTPCLikelihoodMerge.hxx"

// Basic algorithm class for the TPC reconstruction. It takes cares of the input-output data handling and calling the different algorithms in the right order.
//

class TRexReco{
public:
  
  // Default constructor
  TRexReco();
  // Default destructor.
  virtual ~TRexReco();
  
  // Default algorithm process method. This is the core of the tpcRecon algorithm.
  ND::THandle<ND::TAlgorithmResult> Process(const ND::TAlgorithmResult&);
  
  // Set the diffusion Sigma for the likelihood fit.
  // Well, the fitter should do that. This is the property of the fit
  /*
    void SetSigma(double sigma0,double sigma1,double sigma2){
    Sigma0 = sigma0;
    SigmaZ = sigma1;
    SigmaZ2 = sigma2;
    likfitparam->SetSigma(Sigma0,SigmaZ,SigmaZ2); // Comunicate it to the likfit
    }*/
  
private: 
  
  /// Seeding
  ND::TTPCSeeding* fSeedingAlgo;

  /// T0 finder
  ND::TTPCT0Finder* fT0Algo;

  /// Tracking algorithm taking care of likelihood fit
  /// and other things like momentum by range.
  ND::TTPCTracking* fTrackingAlgo;

  /// Calculation of the matching chi2 between tracks
  /// connected by a junction.
  ND::TTPCLikelihoodMatch* fLklhdMatchAlgo;

  /// PID algorithm taking care of PID ... Duh !
  ND::TTPCPID* fPIDAlgo;

  /// Merge paths crossing a junction or from a broken track
  ND::TTPCLikelihoodMerge* fLikelihoodMergeAlgo;

};

#endif
