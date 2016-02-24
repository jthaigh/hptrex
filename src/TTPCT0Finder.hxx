#ifndef TTPCT0Finder_hxx_seen
#define TTPCT0Finder_hxx_seen

#include <TND280Event.hxx>
#include <THit.hxx>
#include <THitSelection.hxx>
#include <TH2F.h>
//#include <HEPUnits.hxx>

#include "TTPCT0.hxx"
#include "TTPCPattern.hxx"

namespace ND {
  class TTPCT0Finder;
}

/// Algorithm to determine the T0 of the TPC tracks based on the matching with the other detectors.
class ND::TTPCT0Finder {
  public:
    /// Default constructor.
    TTPCT0Finder();
    /// Default destructor.
    virtual ~TTPCT0Finder();

    /// Calculate a default T0 using all the hits from all the patterns.
    void CalculateDefaultT0(ND::TReconObjectContainer *allPatterns);

    /// Perform the T0 search
    void Process(ND::THandle<ND::TTPCPattern> pattern);

    /// Prepare hits for fit
    void PrepareScintHits(ND::TND280Event& event);

    /// Find the best T0 for a cathode crossing track.
    /// Either pick the T0 from one of the segment or calculate the T0 from hits at the cathode.
    void FindCathodeCrosserT0(ND::THandle<ND::TTPCPath> PathA, ND::THandle<ND::TTPCPath> PathB, TTPCT0 &T0found);


  private:

    /// DRift velocity from the data base.
    double fDriftVelocity;
    /// ECAL enable flag.
    bool enabledECAL;
    /// FGD enable flag.
    bool enabledFGD;
    /// P0D enable flag.
    bool enabledP0D;
    /// SMRD enable flag.
    bool enabledSMRD;

    /// THitSelection with prepared P0D hits
    ND::THitSelection fP0D;

    /// THitSelection with prepared FGD hits
    ND::THitSelection fFGD_TPC1; //relevant for TPC1
    ND::THitSelection fFGD_TPC2; //relevant for TPC2
    ND::THitSelection fFGD_TPC3; //relevant for TPC3

    /// THitSelection with prepared ECAL hits
    ND::THitSelection fECAL_TPC1; //relevant for TPC1
    ND::THitSelection fECAL_TPC2; //relevant for TPC2
    ND::THitSelection fECAL_TPC3; //relevant for TPC3
    ND::THitSelection fECAL_SIDE; //relevant from side

    /// THitSelection with prepared SMRD hits
    ND::THitSelection fSMRD;

    /// Return the time computed from the seed values (y0,z0,ty0,rho)
    bool FindT0UsingSeed(ND::THandle<ND::TTPCPath> Path, TTPCT0 &T0found);

    /// Return the time computed by looking at the spread in time
    /// of the hits in the pattern
    void FindTimeRange(ND::THandle<ND::TTPCPattern> Pattern, double &minTime, double &maxTime);

    /// Return the t0 computed by looking at the spread in time
    /// of the hits in the pattern
    bool FindT0Range(ND::THandle<ND::TTPCPattern> Pattern, TTPCT0 &T0found);

    /// Enable the usage of FGD as reference.
    void enableFGD() { enabledFGD = true;}
    /// Disable the usage of FGD as reference.
    void disableFGD() { enabledFGD = false;}

    /// Enable the usage of ECAL as reference.
    void enableECAL() { enabledECAL = true;}
    /// Disable the usage of ECAL as reference.
    void disableECAL() { enabledECAL = false;}

    /// Enable the usage of P0D as reference.
    void enableP0D() { enabledP0D = true;}
    /// Disable the usage of P0D as reference.
    void disableP0D() { enabledP0D = false;}

    /// Enable the usage of SMRD as reference.
    void enableSMRD() { enabledSMRD = true;}
    /// Disable the usage of P0D as reference.
    void disableSMRD() { enabledSMRD = false;}

// TODO: Once the cluster containers are available, this method should be rewritten
//    /// Find the t0 for tracks crossing the cathode.
//    bool findT0fromCrossing( std::vector<TPCMeas> &planes,double &t0, bool halfOfCrossingTrk = false);


    /// Algorithm to obtain the best T0 from FGD.
    /// Call this after calling PrepareHits just to study the T0 determination.
    void fitT0fromFGD(ND::THitSelection fgd, State& seedState, TTPCT0 &thisT0);
    /// Algorithm to obtain the best T0 from ECAL.
    /// Call this after calling PrepareHits just to study the T0 determination.
    void fitT0fromECAL(ND::THitSelection ecal, State& seedState, TTPCT0 &thisT0);
    /// Algorithm to obtain the best T0 from P0D.
    /// Call this after calling PrepareHits just to study the T0 determination.
    void fitT0fromP0D(ND::THitSelection fP0D, const State& seedState, TTPCT0 &thisT0);
    /// Algorithm to obtain the best T0 from SMRD.
    /// Call this after calling PrepareHits just to study the T0 determination.
    void fitT0fromSMRD(ND::THitSelection fSMRD, const State& seedState, TTPCT0 &thisT0);

    /// Internal method propagating the seed to the given hits to find the best match
    void matchStateToHit(const ND::THitSelection& hits, const State& sedd, TTPCT0 &T0Res);

    /// Internal method to fill histograms about the T0 results
    void FillT0Histos(int bin, int tpc, TTPCT0 &T0Res);

  public:

    /// Prepare FGD hits to speed up the processing.
    void prepareFGDHits(ND::TND280Event& event);
    /// Prepare P0D hits to speed up the processing.
    void prepareP0DHits(ND::TND280Event& event);
    /// Prepare ECAL hits to speed up the processing.
    void prepareECALHits(ND::TND280Event& event);
    /// Prepare SMRD hits to speed up the processing.
    void prepareSMRDHits(ND::TND280Event& event);
    /// Empty the containers of hits from the scintillator
    /// detectors to avoid memory leaks.
    void CleanHits();


  private:
    // Length of the TPC drift volume
    double fMax_Drift;

    // parameters for hit charge cuts for hits
    // used in t0 determination
    double fFGD_Qmin;
    double fECAL_Qmin;
    double fP0D_Qmin;
    double fSMRD_Qmin;
    double fFGD_YedgeSize;
    int    fFGD_Layers;
    int    fECAL_Layers;
    int    fP0D_Layers;

    /// chi2 cut maximum to be used for t0 determination
    double fChi2_Max;
    /// Prevent crazy matches
    double fSafetyResidualCut;

    // diagnostic histograms
    bool fEnableHistos;
    static TH1F* hFGDt0;
    static TH1F* hECALt0;
    static TH1F* hP0Dt0;
    static TH1F* hSMRDt0;
    static TH1F* hUsedt0;
    static TH1F* hChi2t0;
    static TH1F* hChi2t0nolog;
    static TH1F* hUsedt0From;
    static TH1F* hCathodet0;

    static TH1F* hFGDmP0Dt0;
    static TH1F* hFGDmECALt0;
    static TH1F* hP0DmECALt0;

    static TH2F* hPositiontpc1t0;
    static TH2F* hPositiontpc2t0;
    static TH2F* hPositiontpc3t0;

};


#endif
