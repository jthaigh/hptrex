#include "TTPCTracking.hxx" 

#include "TTPCCalibration.hxx" 
#include "TTPCHVCluster.hxx" 
#include "TTPCUtils.hxx" 

#include <TOARuntimeParameters.hxx>
#include <TND280Event.hxx>
#include <TEventFolder.hxx>

//*****************************************************************************
trex::TTPCTracking::TTPCTracking( ){
//*****************************************************************************

  fMinNumberOfClusters  = 5;

  fRunLikelihoodFit=true;

  fLklhdFitPath = new trex::TTPCLikFitPath();

}

//*****************************************************************************
trex::TTPCTracking::~TTPCTracking( ){
//*****************************************************************************
  delete fLklhdFitPath;
}

//*****************************************************************************
void trex::TTPCTracking::Process(trex::TTRExPattern& Pattern){
//*****************************************************************************
  // return false when no seed was found for any path ???
  std::vector< trex::TTRExPath >& Paths = Pattern.GetPaths();
  
  for (auto pth = Paths.begin(); pth != Paths.end(); pth++) {
    trex::TTRExPath& path = *pth;
    // Having a seed is a requirement to the likelihood fit
    if (fRunLikelihoodFit && path.CheckStatus(trex::kChi2Fit) && !path->CheckStatus(trex::kRan)){ 
      LikelihoodFit( path);
    }

  }
}


// *****************************************************************************
void trex::TTPCTracking::LikelihoodFit(trex::TTRExPath>& thePath){
// *****************************************************************************
  // Clear status bit
  thePath.ClearStatus(trex::kLikelihoodFit);
  thePath.SetStatus(trex::kRan);

  std::vector<TTRExHVCluster>& allClusters = thePath.GetHits();

  if( allClusters->size() == 0 ) {
    return; 
  }

  //MDH TODO: Check if this needs to store clusters by reference.
  std::vector<TTRExHVCluster> selectedClu;
  std::vector<double> seedState = thePath.GetFrontSeedState();
  fLklhdFitPath->Reset();
  fLklhdFitPath->PrepareClustersForFitting(allClusters, selectedClu, seedState[3]);

  fLklhdFitPath->SetupLogLklhdMinimizer(thePath);

  // Returns 0 when Minuit succeeded
  if (fLklhdFitPath->LogLklhdMinimizer(selectedClu) == 0){
    fLklhdFitPath->SaveFitResults(thePath);
  }

  fLklhdFitPath->Reset();

  // We can use the truth as fake seed to check the impact of the seeding

  //MDH TODO: Is this still needed?
  // Important to prevent memory leaks
  if( selectedClu->size() > 0 ) selectedClu->erase(selectedClu->begin(),selectedClu->end());

}
