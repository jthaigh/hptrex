#include "TTPCTracking.hxx" 

#include "TTRExHVCluster.hxx" 
#include "TTPCUtils.hxx" 

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
    std::cout<<"Tracking sees a path, HasChi2Fit="<<path.HasChi2Fit()<<", HasRunFit="<<path.HasRunFit()<<std::endl;
    if (fRunLikelihoodFit && path.HasChi2Fit() && !path.HasRunFit()){ 
      LikelihoodFit( path);
    }

  }
}


// *****************************************************************************
void trex::TTPCTracking::LikelihoodFit(trex::TTRExPath& thePath){
// *****************************************************************************
  // Clear status bit
  thePath.SetHasLikelihoodFit(false);
  thePath.SetHasRunFit(true);

  std::vector<TTRExHVCluster*>& allClusters = thePath.GetClusters();

  if( allClusters.size() == 0 ) {
    return; 
  }

  std::vector<TTRExHVCluster*> selectedClu;
  std::vector<double> seedState = thePath.GetFrontSeedState();
  fLklhdFitPath->Reset();
  fLklhdFitPath->PrepareClustersForFitting(allClusters, selectedClu, seedState[3]);

  fLklhdFitPath->SetupLogLklhdMinimizer(thePath);

  // Returns 0 when Minuit succeeded
  std::cout<<"Performing likelihood fit..."<<std::endl;
  if (fLklhdFitPath->LogLklhdMinimizer(selectedClu) == 0){
    std::cout<<"Saving results..."<<std::endl;
    fLklhdFitPath->SaveFitResults(thePath);
  }
    std::cout<<"Succeeded!"<<std::endl;
  fLklhdFitPath->Reset();

}
