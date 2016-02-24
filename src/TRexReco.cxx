#include <cmath>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <memory>

#include <TND280Event.hxx>
#include <TGeomInfo.hxx>
#include <TOARuntimeParameters.hxx>
#include <THandle.hxx>
#include <TAlgorithmResult.hxx>

#include "TRexReco.hxx"
#include "TTPCCalibration.hxx"
#include "TTPCRecPackUtils.hxx"
#include "TTPCUtils.hxx"


/*******************
 ** Constructor   **
 ******************/
TRexReco::TRexReco(): ND::TAlgorithm("TRexReco") {

  ND::TND280Log::SetLogLevel("TREx",ND::TND280Log::QuietLevel);

  TTPCRecPackUtils::InitRecPackManager();

  // Output
  fSavePatRecResults = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.SavePatternRecognitionResults");
  fSaveUnmergedResults = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.SaveGasInteractionOutput");
  fSaveTRExContainers = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.SaveTRExContainers");

  // skip complicated shower events
  fMaxNbWaveforms = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.MaxNbWaveforms");

  // Waveform treatment
  fWaveformAlgo = new TTPCWaveformDigest();

  // Pattern reco
  fPatRecAlgo = new ND::TTPCTRExPatAlgorithm();

  // Seeding
  fSeedingAlgo = new ND::TTPCSeeding();

  // T0 fit 
  fT0Algo = new ND::TTPCT0Finder();

  // Merge across the MM horizontal gap
  fMMHoriGapMMergeAlgo = new ND::TTPCMMHoriGapMerge();
  // Merge across the MM vertical gap
  fMMVertGapMMergeAlgo = new ND::TTPCMMVertGapMerge();

  // Tracking: Likelihood fit and others like momentum by range
  fTrackingAlgo = new ND::TTPCTracking();

  // Matching likelihood between tracks
  fLklhdMatchAlgo = new ND::TTPCLikelihoodMatch();

  // Pid
  fPIDAlgo = new ND::TTPCPID();

  // Merge all Likelihood matches
  fLikelihoodMergeAlgo = new ND::TTPCLikelihoodMerge();

  // Merge across the MM vertical gap
  fCathCrosserMergeAlgo = new ND::TTPCCathCrosserMerge(fT0Algo);

  //ftimer.start(NULL);

}


/*******************
 ** Destructor   **
 ******************/
TRexReco::~TRexReco() {

  delete fWaveformAlgo;
  delete fPatRecAlgo;
  delete fSeedingAlgo;

  delete fT0Algo;
  delete fMMHoriGapMMergeAlgo;
  delete fMMVertGapMMergeAlgo;
  delete fTrackingAlgo;
  delete fLklhdMatchAlgo;
  delete fPIDAlgo;
  delete fLikelihoodMergeAlgo;
  delete fCathCrosserMergeAlgo;

}




/*******************************************************************
 ** Main function that calls each algorithm and controls em all   **
 *******************************************************************/
ND::THandle<ND::TAlgorithmResult> TRexReco::Process(const ND::TAlgorithmResult& in) {

  ND::THandle<ND::TAlgorithmResult> result( CreateResult() );
  ND::THitSelection* Used = new ND::THitSelection("TRExused");
  ND::THitSelection* Unused = new ND::THitSelection("TRExunused");

  unsigned int stepCount = 1;

  //ftimer.restart(NULL);
  ND::TND280Event& event = GetEvent();


  //___________ Read the calibration tables ___________
  ND::tpcCalibration().ReadCalibration(event);

  
  //___________ Hit preparation ___________
  if ( ND::tpcDebug().GeneralSteps(DB_INFO))
    std::cout<<" ============= TRExReco step "<<stepCount++<<": Hit preparation"<<std::endl;
  ND::THandle<ND::THitSelection> tpc(new ND::THitSelection("tpc"));
  fWaveformAlgo->Process(event, tpc);

  if( tpc->size() > fMaxNbWaveforms){
    std::cout<< "TREX WARNING: too many waveforms in the TPCs ("<<tpc->size()<<" waveforms). Give up the TPC reconstruction entirely"<<std::endl;
    return result;
  }
  
  //___________ Pattern recognition ___________
  // These patterns are for each MM column of each read out plane independently
  // Save the pattern recognition results if we want to.
  if ( ND::tpcDebug().GeneralSteps(DB_INFO))
    std::cout<<" ============= TRExReco step "<<stepCount++<<": Pattern recognition"<<std::endl;
  fPatRecAlgo->Process(tpc, Used, Unused);
  result->AddHitSelection(Used);
  result->AddHitSelection(Unused);
  ND::TReconObjectContainer *PartColumnPatterns = new ND::TReconObjectContainer("TRExPatternReco");
  fPatRecAlgo->GetPatterns(PartColumnPatterns);

  // Define a T0 that puts all the hits in the drift volume if possible.
  fT0Algo->CalculateDefaultT0(PartColumnPatterns);

  if ( ND::tpcDebug().GeneralSteps(DB_VERBOSE)) {
    PrintoutObjects(PartColumnPatterns);
  }

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


  //___________ Merging across MM horizontal gap ___________
  // Some of the patterns are broken at the horizontal gap.
  // If a pattern is not matched to another one, then it is just copied to ColumnPatterns.
  // If two pattern are matched, only one is saved to the output with the added content of the second pattern.
  // The second pattern can be deleted ... we have to be careful in terms of memory management.
  if ( ND::tpcDebug().GeneralSteps(DB_INFO))
    std::cout<<" ============= TRExReco step "<<stepCount++<<": MM horizontal gap merging"<<std::endl;
  ND::TReconObjectContainer *ColumnPatterns = new ND::TReconObjectContainer("TRExPatternReco");
  fMMHoriGapMMergeAlgo->Process(PartColumnPatterns, ColumnPatterns);

  if ( ND::tpcDebug().GeneralSteps(DB_VERBOSE)) {
    PrintoutObjects(ColumnPatterns);
  }

  //___________ Merging across MM vertical gap ___________
  // Now our patterns are for each read out plane.
  // If a pattern is not matched to another one, then it is just copied to ReadOutPatterns.
  // If two pattern are matched, only one is saved to the output with the added content of the second pattern.
  // The second pattern can be deleted ... we have to be careful in terms of memory management.
  if ( ND::tpcDebug().GeneralSteps(DB_INFO))
    std::cout<<" ============= TRExReco step "<<stepCount++<<": MM vertical gap merging"<<std::endl;
  ND::TReconObjectContainer *ReadOutPatterns = new ND::TReconObjectContainer("TRExUnmergedReco");
  fMMVertGapMMergeAlgo->Process(ColumnPatterns, ReadOutPatterns);

  if ( ND::tpcDebug().GeneralSteps(DB_VERBOSE)) {
    PrintoutObjects(ReadOutPatterns);
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

  //___________ Merging tracks crossing the TPC ___________
  if ( ND::tpcDebug().GeneralSteps(DB_INFO))
  std::cout<<" ============= TRExReco step "<<stepCount++<<": Merging through going and broken tracks"<<std::endl;
  ND::TReconObjectContainer *ReadOutMergedPatterns = new ND::TReconObjectContainer("TRExLklhdMergedReco");
  fLikelihoodMergeAlgo->Process(ReadOutPatterns, ReadOutMergedPatterns);

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

  if ( ND::tpcDebug().GeneralSteps(DB_VERBOSE)) {
    PrintoutObjects(ReadOutMergedPatterns);
  }

  //___________ Merging across cathode ___________
  if ( ND::tpcDebug().GeneralSteps(DB_INFO))
  std::cout<<" ============= TRExReco step "<<stepCount++<<": Merging cathode crossers"<<std::endl;
  ND::TReconObjectContainer *FullTPCPatterns = new ND::TReconObjectContainer("TRExAllMergedReco");
  fCathCrosserMergeAlgo->Process(ReadOutMergedPatterns, FullTPCPatterns);

  if ( ND::tpcDebug().GeneralSteps(DB_VERBOSE)) {
    PrintoutObjects(FullTPCPatterns);
  }

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


// *******************************************************************
void TRexReco::StoreContainerInOutput( ND::THandle<ND::TAlgorithmResult> algoResult, ND::TReconObjectContainer *container, std::string oaEventName){
  if (fSaveTRExContainers){
    algoResult->AddResultsContainer(container);
  }
  else {
    // Convert TTPCPatterns to TReconVertex if there is one or more Junction.
    // Convert to TReconTrack if there is no Junction.
    // This is done automatically.
    ND::TReconObjectContainer *OutputPatterns = new ND::TReconObjectContainer(oaEventName.c_str());
    for (ND::TReconObjectContainer::iterator pattern = container->begin(); pattern != container->end(); pattern++) {
      ND::THandle<ND::TTPCPattern> Pattern = (*pattern);
      ND::THandle<ND::TReconBase> oaEvtPattern = Pattern->ConvertToOAEvent();
      if ( ND::tpcDebug().GeneralSteps(DB_VVERBOSE)) {
        std::cout<<" ============= Constituent structure"<<std::endl;
        TTPCUtils::PrintConstituentMap(oaEvtPattern,0);
      }
      OutputPatterns->push_back(oaEvtPattern);
    }
    algoResult->AddResultsContainer(OutputPatterns);
  }

}


// *******************************************************************
void TRexReco::PrintoutObjects(ND::TReconObjectContainer *inPatterns){
  for (ND::TReconObjectContainer::iterator colPatIt = inPatterns->begin(); colPatIt < inPatterns->end(); colPatIt++)
  {
    ND::THandle<ND::TTPCPattern> pattern = *colPatIt;
    if (!pattern){
      std::cout << "   ============= TRExReco SAFETY WARNING: Object from pattern recognition is NOT a TTPCPattern !" << std::endl;
      continue;
    }

    ND::THandle<ND::TReconObjectContainer> constituents = pattern->GetConstituents();
    int nconst = constituents->size();
    if (!nconst) {
      std::cout << " ============= TRExReco WARNING: TTPCPattern has 0 constituents !" << std::endl;
      continue;
    }
    std::cout << " =============> pattern, id "<<pattern->GetId() << std::endl;

    for (ND::TReconObjectContainer::iterator constIt = constituents->begin(); constIt != constituents->end(); constIt++)
    {
      ND::THandle<ND::TTPCPath> path = *constIt;
      ND::THandle<ND::TTPCJunction> junction = *constIt;
      if (path){
        std::cout << "  |-> path, id "<<path->GetId() << std::endl;
      } else if (junction) {
        std::cout << "  |-> junction, id "<<junction->GetId() << std::endl;
        for (ND::TReconObjectContainer::iterator pathIt = junction->GetConstituents()->begin(); pathIt != junction->GetConstituents()->end(); pathIt++){
          path = *pathIt;
          std::cout << "  |   |-> path, id "<<path->GetId() << std::endl;
        }
      } else {
        std::cout << " ============= TRExReco WARNING: A constituent of the pattern is neither a TTPCPath, nor a TTPCJunction !" << std::endl;
      }
    }
  }
}
