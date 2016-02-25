#include "TTPCTracking.hxx" 

#include "TTPCCalibration.hxx" 
#include "TTPCHVCluster.hxx" 
#include "TTPCUtils.hxx" 

#include <TOARuntimeParameters.hxx>
#include <TND280Event.hxx>
#include <TEventFolder.hxx>

//*****************************************************************************
ND::TTPCTracking::TTPCTracking( ){
//*****************************************************************************

  fRunLikelihoodFit = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.RunLikelihoodFit");

  fUseTruthAsFitResult = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.UseTruthAsFitResults");

  fMinNumberOfClusters  = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.LikFit.MinNbClusters");

  fLklhdFitPath = new TTPCLikFitPath();

  fClusterCorrection = new TTPCClusterCorrection();
}

//*****************************************************************************
ND::TTPCTracking::~TTPCTracking( ){
//*****************************************************************************
  delete fLklhdFitPath;
}

//*****************************************************************************
void ND::TTPCTracking::Process(ND::THandle<ND::TTPCPattern> Pattern){
//*****************************************************************************
  // return false when no seed was found for any path ???
  std::vector< ND::THandle<ND::TTPCPath> > Paths = Pattern->GetPaths();
  
  for (std::vector< ND::THandle<ND::TTPCPath> >::iterator pth = Paths.begin(); pth != Paths.end(); pth++) {
    ND::THandle<ND::TTPCPath> path = *pth;
    // Having a seed is a requirement to the likelihood fit
    if (fRunLikelihoodFit && path->CheckStatus(ND::TReconBase::kChi2Fit) && !path->CheckStatus(ND::TReconBase::kRan)){ 
      LikelihoodFit( path);
    }

  }
}


// *****************************************************************************
void ND::TTPCTracking::LikelihoodFit(ND::THandle<ND::TTPCPath> thePath){
// *****************************************************************************
  // Clear status bit
  thePath->ClearStatus(ND::TReconBase::kLikelihoodFit);
  thePath->SetStatus(ND::TReconBase::kRan);

  ND::THandle<ND::THitSelection> allClusters = thePath->GetHits();
  if( !allClusters ) {
    if( ND::tpcDebug().Tracking(DB_INFO) )
      std::cout << " TTPCLikFitPath: NO HITS in track, abort fitting. " << std::endl; 
    return; 
  }
  if( allClusters->size() == 0 ) {
    if( ND::tpcDebug().Tracking(DB_INFO) )
      std::cout << " TTPCLikFitPath: NO HITS in track, abort fitting. " << std::endl; 
    return; 
  }

  // Apply various field/empirical corrections to the clusters
  fClusterCorrection->Apply(thePath);

  ND::THandle<ND::THitSelection> selectedClu(new THitSelection());
  State seedState = thePath->GetFrontSeedState();
  fLklhdFitPath->Reset();
  fLklhdFitPath->PrepareClustersForFitting(allClusters, selectedClu, seedState.vector()[3]);



  fLklhdFitPath->SetupLogLklhdMinimizer(thePath);

  // Returns 0 when Minuit succeeded
  if (fLklhdFitPath->LogLklhdMinimizer(selectedClu) == 0){
    fLklhdFitPath->SaveFitResults(thePath);
  }
  // TODO: Refits with different E or B field corrections
  /*
  for ( unsigned int cf = 1; cf < fClusterCorrection->GetNbConfigurations(); cf++){
    fLklhdFitPath->ResetMinuitParam();
    fClusterCorrection->LoadConfiguration(cf);
    if (fLklhdFitPath->Process(selectedClu, Result) == 0){
      thePath->SaveAlternateFitState(Result);
    }
  }
  */

  fLklhdFitPath->Reset();


  if ( ND::tpcDebug().Tracking(DB_INFO)){
    std::cout<<" =========================================================================="<<std::endl;
    std::cout<<"  Path Id "<<thePath->GetId()<<std::endl;
    if( thePath->HasFitState()){
      State fitResult = thePath->GetFitState();
      std::cout<<"  Likelihood: "<<fitResult.vector()<<std::endl;
    } else {
      std::cout<<"  Likelihood fit failed !"<<std::endl;
    }
  }

  // We can use the truth as fake seed to check the impact of the seeding
  if (ND::tpcCalibration().IsMC() && 
    ( fUseTruthAsFitResult || ND::tpcDebug().Tracking(DB_INFO) )) {
      State trueState = State(7);
      // This gets the truth or returns false. Includes check if input file is MC.
      bool ok = TTPCUtils::TrueStateNear3Dpoint(thePath->GetHits(), thePath->GetFirstPosition(), trueState);
      if (ok && fUseTruthAsFitResult)
        thePath->SaveFitState(trueState);
      if (ok && ND::tpcDebug().Tracking(DB_INFO))
        std::cout<<"  Truth: "<<trueState.vector()<<std::endl;
  }

  if ( ND::tpcDebug().Tracking(DB_INFO)){
    std::cout<<" ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  }

  // Important to prevent memory leaks
  if( selectedClu->size() > 0 ) selectedClu->erase(selectedClu->begin(),selectedClu->end());

}
