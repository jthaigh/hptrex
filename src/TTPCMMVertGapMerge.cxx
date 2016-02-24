#include "TTPCMMVertGapMerge.hxx"
#include "TTPCJunction.hxx"
#include "TTPCRecPackUtils.hxx"
#include "TTPCCalibration.hxx"
#include "TTPCDebug.hxx"

#include <THandle.hxx>
#include <TGeomInfo.hxx>
#include <TOARuntimeParameters.hxx>



//*****************************************************************************
ND::TTPCMMVertGapMerge::TTPCMMVertGapMerge(void) : ND::TTPCSeed2CluMatchMerge(){
  fRun = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.RunMMVerticalGapMerging");
  if (!fRun)
    std::cout<<"TRexRecon WARNING: MM vertical gap matching and merging disabled"<<std::endl;

  fMinDistForMatchP2J = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.MMVertGapMerge.MinDistForMatchPathToJunction");
  fMinDistForMatchP2P = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.MMVertGapMerge.MinDistForMatchPathToPath");
  fNbWantedColumns = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.MMVertGapMerge.NbColumnsToGap");

  fAlgo = kMMVertGapMerge;

  fFakeJunctionForStudies = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.MMVertGapMerge.FakeJunctionForStudies");
}



//*****************************************************************************
void ND::TTPCMMVertGapMerge::Process(ND::TReconObjectContainer *inputPatterns, ND::TReconObjectContainer *mergedPatterns){
  if ( !fRun){
    for (ND::TReconObjectContainer::iterator pattit = inputPatterns->begin(); pattit != inputPatterns->end(); pattit++) {
      mergedPatterns->push_back(*pattit);
    }
    return;
  }

  if ( ND::tpcDebug().MMVertGapMerge(DB_INFO))
    std::cout<<"TTPCMMVertGapMerge::Process"<<std::endl;

  std::vector< ND::THandle<ND::TReconBase> > rawPatterns[3][2][2];
  for (ND::TReconObjectContainer::iterator pattit = inputPatterns->begin(); pattit != inputPatterns->end(); pattit++) {
    ND::THandle<ND::TTPCPattern> pattern = *pattit;
    int tpc = -1;
    int rp = -1;
    int updown = -1;
    FindLocation(pattern, tpc, rp, updown);
    pattern->SetUsable(true);
    rawPatterns[tpc][rp][updown].push_back(pattern);
    if ( ND::tpcDebug().MMVertGapMerge(DB_VERBOSE)){
      std::cout<<" Input pattern: TPC"<<(tpc+1)<<", RP"<<rp<<", is downstream = "<<updown<<std::endl;
      for (ND::TReconObjectContainer::iterator constit = pattern->GetConstituents()->begin(); constit != pattern->GetConstituents()->end(); constit++) {
        ND::THandle<ND::TTPCPath> tmpPath = (*constit);
        if (!tmpPath)
          continue;
        std::cout<<"        Path Id "<<tmpPath->GetId()<<std::endl;
        ND::THandle<ND::TTPCHVCluster> tmpClu = *(tmpPath->GetHits()->begin());
        std::cout<<"        First cluster at "<<tmpClu->GetPosition().X()<<"   "<<tmpClu->GetPosition().Y()<<"   "<<tmpClu->GetPosition().Z()<<std::endl;
        tmpClu = *(tmpPath->GetHits()->rbegin());
        std::cout<<"        Last cluster at "<<tmpClu->GetPosition().X()<<"   "<<tmpClu->GetPosition().Y()<<"   "<<tmpClu->GetPosition().Z()<<std::endl;
      }
    }
  }

  fMaxNbMatchCand = 0;
  fMaxNbPropagCand = 0;
  // Loop through the TPCs and Readout planes.
  //   Loop through the patterns and find the candidates for matching and propagating.

  for( int tpc = 0; tpc < 3; tpc++){
    for( int rp = 0; rp < 2; rp++){
      if ( ND::tpcDebug().MMVertGapMerge(DB_VERBOSE))
        std::cout<<" ==> Match Merge for TPC "<<(tpc+1)<<", RP "<<rp<<std::endl;
      int UpIdx = 0;
      int DownIdx = 1;
      std::vector< ND::THandle<ND::TReconBase> > resultPatterns;

      // ---------- Propagate the downstream paths to upstream patterns
      // Clear the containers
      CleanUp();
      // Fill the containers for the downstream paths that I will propagate to upstream patterns
      for (ND::TReconObjectContainer::iterator pattit = rawPatterns[tpc][rp][DownIdx].begin(); pattit != rawPatterns[tpc][rp][DownIdx].end(); pattit++) {
        ND::THandle<ND::TTPCPattern> pattern = *pattit;
        FindPropagCandidates(pattern, DownIdx);
      }

      // Process the upstream patterns
      if ( ND::tpcDebug().MMVertGapMerge(DB_VERBOSE))
        std::cout<<"   Process the upstream patterns"<<std::endl;
      MatchPatternsAWithB(rawPatterns[tpc][rp][UpIdx], UpIdx, DownIdx, resultPatterns);

      // ---------- Propagate the upstream paths to downstream patterns
      // Clear the containers
      fNbMatchCand = fMaxNbMatchCand;
      fNbPropagCand = fMaxNbPropagCand;
      CleanUp();
      // Fill the containers for the upstream paths that I will propagate to downstream patterns
      for (ND::TReconObjectContainer::iterator pattit = resultPatterns.begin(); pattit != resultPatterns.end(); pattit++) {
        ND::THandle<ND::TTPCPattern> pattern = *pattit;
        pattern->SetUsable(true);
        FindPropagCandidates(pattern, UpIdx);
      }

      // Process the downstream patterns
      if ( ND::tpcDebug().MMVertGapMerge(DB_VERBOSE))
        std::cout<<"   Process the downstream patterns"<<std::endl;
      MatchPatternsAWithB(rawPatterns[tpc][rp][DownIdx], DownIdx, UpIdx, resultPatterns);

      for (ND::TReconObjectContainer::iterator pattit = resultPatterns.begin(); pattit != resultPatterns.end(); pattit++) {
        ND::THandle<ND::TTPCPattern> pattern = *pattit;
        if ( ! pattern->IsUsable())
          continue;
        mergedPatterns->push_back(pattern);
      }
    }
  }
  
  // Clean up to avoid memory leaks and so on.
  fNbMatchCand = fMaxNbMatchCand;
  fNbPropagCand = fMaxNbPropagCand;
  CleanUp();
}


//*****************************************************************************
void ND::TTPCMMVertGapMerge::FindLocation(ND::THandle<ND::TTPCPattern> pattern, int &tpc, int &rp, int &updown){
  ND::TReconObjectContainer::iterator constit = pattern->GetConstituents()->begin();
  ND::THandle<ND::TTPCHitPad> hitPad = *((*constit)->GetHits()->begin());
  ND::THandle<ND::TTPCHVCluster> Cluster = *((*constit)->GetHits()->begin()); 
  if ( Cluster) {
    hitPad = *(Cluster->GetHits().begin());
  }
  if (!hitPad){
    // TODO: Proper exception
    throw;
  }
  int geotpc, geohalf, geomm, geopad;
  ND::TTPCGeom& tpcGeom = ND::TGeomInfo::TPC();
  tpcGeom.GetGeometryInfo(hitPad->GetGeomId(), geotpc, geohalf, geomm, geopad);
  tpc = geotpc;
  rp = geohalf;

  int column = ND::TGeomInfo::TPC().ActivePlane(hitPad->GetGeomId());
  // upstream
  if(( -1   < column && column < 36) || 
     ( 71  < column && column <108) ||
     ( 143 < column && column <180) ){
    updown = 0;
  } else if(( 35  < column && column < 72) || 
     (107 < column && column < 144) ||
     (179 < column && column < 216) ){
    updown = 1;
  } else {
    // TODO: Proper exception
    throw;
  }
}





//*****************************************************************************
bool ND::TTPCMMVertGapMerge::HitSelected(ND::THandle<ND::TTPCHitPad> hitPad, int &Zone){
  int column = ND::TGeomInfo::TPC().ActivePlane(hitPad->GetGeomId());
  if(( (35 - fNbWantedColumns) < column && column < 36) || 
     ( (107- fNbWantedColumns) < column && column <108) ||
     ( (179- fNbWantedColumns) < column && column <180) ){
    if ( ! Zone )
      return true;
  } else if(( 35 < column && column < (36 + fNbWantedColumns)) ||
     (107 < column && column < (108 + fNbWantedColumns)) ||
     (179 < column && column < (180 + fNbWantedColumns)) ){
    if ( Zone )
      return true;
  }
  return false;
}
