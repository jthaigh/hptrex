#ifndef TTPCSeeding_hxx_seen
#define TTPCSeeding_hxx_seen

#include "TTRExPattern.hxx" 
#include "TTRExHVCluster.hxx"

#define RHOMIN 5.e-10  // Equivalent to 20 GeV.


namespace trex {
  class TTPCSeeding;
}

class trex::TTPCSeeding {

  public:
  TTPCSeeding();
  ~TTPCSeeding() {};

  void Process(TTRExPattern& Pattern);
    /// This is used by Process but is also called directly
    /// after merging two tracks broken at MM gaps.
  void FindSeed(TTRExPath& thePath);

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

  void PrepareSeeding( std::vector<trex::TTRExHVCluster*>& HVclu );

  void PrepareClustersForSeeding( std::vector<trex::TTRExHVCluster*>& HVclu );
    /// Riemann fit method.
  double Riemann( std::vector<trex::TTRExHVCluster*>& HVclu, std::vector<double>& helixParam);
    /// 3 point fit method.
  double R2( std::vector<trex::TTRExHVCluster*>& HVclu, std::vector<double>& helixParam);

  double CalculateRhoSign(double Y0, double Z0);

  bool fUseTruthAsSeedResult;

  double FinalizeSeed( std::vector<trex::TTRExHVCluster*>& HVclu, double rho, std::vector<double>& finalState);

  bool IsResultValid( std::vector<trex::TTRExHVCluster*>& HVclu, std::vector<double>& Result);
};


#endif 
