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
  for (ND::TReconObjectContainer::iterator pattern = PartColumnPatterns->begin(); pattern != PartColumnPatterns->end(); pattern++) {
      fSeedingAlgo->Process(*pattern);
  }
  for (ND::TReconObjectContainer::iterator pattern = ReadOutPatterns->begin(); pattern != ReadOutPatterns->end(); pattern++) {
    //___________ First tracking pass ___________
    // I didn't call this fitting because I'm thinking that we might want to put here also the momentum by range ... maybe not. Food for thought !
    fTrackingAlgo->Process(*pattern);

    //___________ PID for the reconstructed tracks ___________
    fPIDAlgo->Process(*pattern);

  }
  stepCount+=4;

  //___________ Likelihood calculation between tracks connected by a junction and broken tracks ___________
  fLklhdMatchAlgo->Process(ReadOutPatterns);


  // Here we save the output that will be used in the Prod7 TPC gas interaction analysis.
  if (fSaveUnmergedResults){
    StoreContainerInOutput(result, ReadOutPatterns, "GasInteractionOutput");
  }

  PrintoutObjects(ReadOutPatterns);

  // We can here loop over all the patterns again
  for (ND::TReconObjectContainer::iterator pattit = ReadOutMergedPatterns->begin(); pattit != ReadOutMergedPatterns->end(); pattit++) {
    // The output will be like what is done in tpcRecon in terms of recombining tracks broken by delta rays.
    ND::THandle<ND::TTPCPattern> Pattern = *pattit;
    //___________ Second seeding pass ___________
    fSeedingAlgo->Process(Pattern);
    //___________ Second tracking pass ___________
    // Here we rerun only on tracks that have been merged to create through going tracks.
    // So most tracks will not change here.
    fTrackingAlgo->Process(Pattern);

    //___________ PID for the reconstructed tracks ___________
    fPIDAlgo->Process(Pattern);
    
  }
  stepCount+=3;

  //___________ Save the tpcRecon-like output ___________
  // Here we save the output that will be used in the Prod7 globalRecon.
  StoreContainerInOutput(result, FullTPCPatterns, "TPCPatterns");
  
  // Clean up the RecPack geometry from the propagation surfaces
  // at each cluster.
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


