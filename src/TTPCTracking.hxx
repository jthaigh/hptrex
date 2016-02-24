#ifndef TTPCTracking_hxx_seen
#define TTPCTracking_hxx_seen

#include "TTPCPattern.hxx" 
#include "TTPCPath.hxx"
#include "TTPCLikFitPath.hxx"
#include "TTPCClusterCorrection.hxx"

#define RHOMIN 5.e-10  // Equivalent to 20 GeV.


// Make a processing class: Create one instance of the class instead of one per track.
//   The constructor just reads the parameters needed
//   FindSeed
//       PrepareSeeding takes the TTPCHVClusters of the track and prepares useful quantities
//       Test a vector of seeding algorithm, and chooses the one with the smallest chi2
//       Need structures (HelixModel ?) to store the YZ results and X results
//

namespace ND {
  class TTPCTracking;
}

class ND::TTPCTracking {

  public:
    TTPCTracking();
    ~TTPCTracking();

    void Process(ND::THandle<ND::TTPCPattern> Pattern);

  private:
    void LikelihoodFit(ND::THandle<ND::TTPCPath> thePath);

    TTPCLikFitPath *fLklhdFitPath;
    TTPCClusterCorrection *fClusterCorrection;

    bool fRunLikelihoodFit;
    bool fUseTruthAsFitResult;

	  unsigned int fMinNumberOfClusters;
};


#endif 
