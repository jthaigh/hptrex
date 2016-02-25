#ifndef TTPCSeeding_hxx_seen
#define TTPCSeeding_hxx_seen

#include "TTPCPattern.hxx" 
#include "TTPCPath.hxx"

#define RHOMIN 5.e-10  // Equivalent to 20 GeV.


namespace ND {
  class TTPCSeeding;
}

class ND::TTPCSeeding {

  public:
    TTPCSeeding();
    ~TTPCSeeding() {};

    void Process(ND::THandle<ND::TTPCPattern> Pattern);
    /// This is used by Process but is also called directly
    /// after merging two tracks broken at MM gaps.
    void FindSeed(ND::THandle<ND::TTPCPath> thePath);

  private:
    unsigned int fNbOrientChange;

    double fXfirst;
    double fYfirst;
    double fZfirst;
    double fXlast;
    double fYlast;
    double fZlast; 

    double fYmid;
    double fZmid;

    double fDriftVelocity;  
    double fChi2max;

    /// Exclude clusters with waveforms containing many peaks from the fit.
    bool fExcludeClusterWithManyPeaks;
    /// Exclude clusters containing saturated waveforms from the fit.
    bool fExcludeSaturatedClusters;

    void PrepareSeeding( ND::THandle<ND::THitSelection> HVclu );
    void PrepareClustersForSeeding( ND::THandle<ND::THitSelection> HVclu );

    /// Riemann fit method.
    double Riemann( ND::THandle<ND::THitSelection> HVclu, State &Helix);
    /// 3 point fit method.
    double R2( ND::THandle<ND::THitSelection> HVclu, State &Helix);

    double CalculateRhoSign(double Y0, double Z0);

    bool fUseTruthAsSeedResult;

    double FinalizeSeed( ND::THandle<ND::THitSelection> HVclu, double rho, State &finalState);

    bool IsResultValid( ND::THandle<ND::THitSelection> HVclu, State &Result);
};


#endif 
