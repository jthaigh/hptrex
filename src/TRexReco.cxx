#include <cmath>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <memory>

#include "TRexReco.hxx"
#include "TTPCRecPackUtils.hxx"
#include "TTPCUtils.hxx"


/*******************
 ** Constructor   **
 ******************/
TRexReco::TRexReco(){

  fSeedingAlgo = new ND::TTPCSeeding();

  // T0 fit 
  fT0Algo = new ND::TTPCT0Finder();

  // Tracking: Likelihood fit and others like momentum by range
  fTrackingAlgo = new ND::TTPCTracking();

  // Matching likelihood between tracks
  fLklhdMatchAlgo = new ND::TTPCLikelihoodMatch();

  // Pid
  fPIDAlgo = new ND::TTPCPID();

  // Merge all Likelihood matches
  fLikelihoodMergeAlgo = new ND::TTPCLikelihoodMerge();

}


/*******************
 ** Destructor   **
 ******************/
TRexReco::~TRexReco() {


}




/*******************************************************************
 ** Main function that calls each algorithm and controls em all   **
 *******************************************************************/
ND::THandle<ND::TAlgorithmResult> TRexReco::Process(const ND::TAlgorithmResult& in) {

  fPatRecAlgo->GetPatterns(PartColumnPatterns);

  // Define a T0 that puts all the hits in the drift volume if possible.
  fT0Algo->CalculateDefaultT0(PartColumnPatterns);

  // The pattern recognition didn't find anything. We can stop here.
  if ( ! PartColumnPatterns->size())
  {
    fPatRecAlgo->CleanUp();
    delete PartColumnPatterns;
    return result;
  }


  if (fSavePatRecResults){
    StoreContainerInOutput(result, PartColumnPatterns, "TPCPatternReco");
  }

  //___________ First seeding pass ___________
  // Loop over the patterns and the seeding will take care of looping over the paths.
  // The results of the seeding will be stored in the paths or whatever new container we use. I use TTPCPattern as an example for now.
  if ( ND::tpcDebug().GeneralSteps(DB_INFO))
    std::cout<<" ============= TRExReco step "<<stepCount++<<": First seeding pass"<<std::endl;
  for (ND::TReconObjectContainer::iterator pattern = PartColumnPatterns->begin(); pattern != PartColumnPatterns->end(); pattern++) {
      fSeedingAlgo->Process(*pattern);
  }
  // We don't want to prepare the hits of the other detectors for the T0 search
  // for every single pattern so let's do it here.
  fT0Algo->PrepareScintHits(event);
  // We can here loop over all the patterns in each read out plane
  if ( ND::tpcDebug().GeneralSteps(DB_INFO))
    std::cout<<" ============= TRExReco loop over the patterns for multiple steps on each one"<<std::endl;
  for (ND::TReconObjectContainer::iterator pattern = ReadOutPatterns->begin(); pattern != ReadOutPatterns->end(); pattern++) {
    //___________ T0 determination ___________
    // Pass the entire pattern to the T0 because if the T0 is found for one track, it's valid for the whole pattern.
    // If no T0 is found, we can use some tricks to guess a T0 interval and let the analyzer decide what to do with it.
    if ( ND::tpcDebug().GeneralSteps(DB_INFO))
      std::cout<<" ============= TRExReco step "<<stepCount+1<<": Search for the track T0"<<std::endl;
    fT0Algo->Process(*pattern);

    //___________ First tracking pass ___________
    // I didn't call this fitting because I'm thinking that we might want to put here also the momentum by range ... maybe not. Food for thought !
    if ( ND::tpcDebug().GeneralSteps(DB_INFO))
      std::cout<<" ============= TRExReco step "<<stepCount+2<<": Tracking (likelihood fit)"<<std::endl;
    fTrackingAlgo->Process(*pattern);

    //___________ PID for the reconstructed tracks ___________
    if ( ND::tpcDebug().GeneralSteps(DB_INFO))
      std::cout<<" ============= TRExReco step "<<stepCount+3<<": PID"<<std::endl;
    fPIDAlgo->Process(*pattern);

  }
  stepCount+=4;

  //___________ Likelihood calculation between tracks connected by a junction and broken tracks ___________
  if ( ND::tpcDebug().GeneralSteps(DB_INFO))
    std::cout<<" ============= TRExReco step "<<stepCount++<<": Matching likelihood calculations"<<std::endl;
  fLklhdMatchAlgo->Process(ReadOutPatterns);


  // Here we save the output that will be used in the Prod7 TPC gas interaction analysis.
  if (fSaveUnmergedResults){
    if ( ND::tpcDebug().GeneralSteps(DB_INFO))
      std::cout<<" ============= TRExReco step "<<stepCount++<<": Save results in GasInteractionOutput container"<<std::endl;
    StoreContainerInOutput(result, ReadOutPatterns, "GasInteractionOutput");
  }

  if ( ND::tpcDebug().GeneralSteps(DB_VERBOSE)) {
    PrintoutObjects(ReadOutPatterns);
  }

  // Apply the offset to the pattern/path/junction ids to separate
  // gas and standard output objects.
  ND::tpcCalibration().OffsetIds();

  // We can here loop over all the patterns again
  for (ND::TReconObjectContainer::iterator pattit = ReadOutMergedPatterns->begin(); pattit != ReadOutMergedPatterns->end(); pattit++) {
    // The output will be like what is done in tpcRecon in terms of recombining tracks broken by delta rays.
    ND::THandle<ND::TTPCPattern> Pattern = *pattit;
    //___________ Second seeding pass ___________
    if ( ND::tpcDebug().GeneralSteps(DB_INFO))
    std::cout<<" ============= TRExReco step "<<stepCount<<": Second seeding pass"<<std::endl;
    fSeedingAlgo->Process(Pattern);
    //___________ Second tracking pass ___________
    // Here we rerun only on tracks that have been merged to create through going tracks.
    // So most tracks will not change here.
    if ( ND::tpcDebug().GeneralSteps(DB_INFO))
    std::cout<<" ============= TRExReco step "<<stepCount+1<<": Second tracking (likelihood fit)"<<std::endl;
    fTrackingAlgo->Process(Pattern);

    //___________ PID for the reconstructed tracks ___________
    if ( ND::tpcDebug().GeneralSteps(DB_INFO))
      std::cout<<" ============= TRExReco step "<<stepCount+2<<": PID"<<std::endl;
    fPIDAlgo->Process(Pattern);
    
  }
  stepCount+=3;

  //___________ Save the tpcRecon-like output ___________
  // Here we save the output that will be used in the Prod7 globalRecon.
  if ( ND::tpcDebug().GeneralSteps(DB_INFO))
    std::cout<<" ============= TRExReco step "<<stepCount++<<": Save results in the standard out container"<<std::endl;
  StoreContainerInOutput(result, FullTPCPatterns, "TPCPatterns");
  

  // Clean up the RecPack geometry from the propagation surfaces
  // at each cluster.
  if ( ND::tpcDebug().GeneralSteps(DB_INFO))
    std::cout<<" ============= TRExReco step "<<stepCount++<<": Clean up RecPack surfaces"<<std::endl;
  for (ND::TReconObjectContainer::iterator pattern = ReadOutMergedPatterns->begin(); pattern != ReadOutMergedPatterns->end(); pattern++) {
    ND::THandle<ND::TTPCPattern> Pattern = *pattern;
    std::vector< ND::THandle<ND::TTPCPath> > Paths = Pattern->GetPaths();
    for (std::vector< ND::THandle<ND::TTPCPath> >::iterator pth = Paths.begin(); pth != Paths.end(); pth++) {
      ND::THandle<ND::TTPCPath> Path = (*pth);
      Path->CleanUp();
    }
  }

  //_______Clean up containers to avoid issues with memory ___________
  fPatRecAlgo->CleanUp();
  fT0Algo->CleanHits();

  if (!fSaveTRExContainers)
  {
    delete PartColumnPatterns;
    delete ColumnPatterns;
    delete ReadOutPatterns;
    delete ReadOutMergedPatterns;
    delete FullTPCPatterns;
  }

  return result;
}


