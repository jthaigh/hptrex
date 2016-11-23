#ifndef TTPCTracking_hxx_seen
#define TTPCTracking_hxx_seen

#include "TTRExPattern.hxx" 
#include "TTRExPath.hxx"
#include "TTPCLikFitPath.hxx"

#define RHOMIN 5.e-10  // Equivalent to 20 GeV.


// Make a processing class: Create one instance of the class instead of one per track.
//   The constructor just reads the parameters needed
//   FindSeed
//       PrepareSeeding takes the TTPCHVClusters of the track and prepares useful quantities
//       Test a vector of seeding algorithm, and chooses the one with the smallest chi2
//       Need structures (HelixModel ?) to store the YZ results and X results
//

namespace trex {
  class TTPCTracking;
}

class trex::TTPCTracking {

  public:
    TTPCTracking();
    ~TTPCTracking();

    void Process(trex::TTRExPattern& Pattern);

  private:
    void LikelihoodFit(trex::TTRExPath& thePath);

    TTPCLikFitPath *fLklhdFitPath;

    bool fRunLikelihoodFit;

    unsigned int fMinNumberOfClusters;
};


#endif 
