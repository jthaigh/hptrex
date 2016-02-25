#include "TTPCMMHoriGapMerge.hxx"
#include "TTPCJunction.hxx"
#include "TTPCRecPackUtils.hxx"
#include "TTPCCalibration.hxx"
#include "TTPCDebug.hxx"

#include <THandle.hxx>
#include <TGeomInfo.hxx>
#include <TOARuntimeParameters.hxx>



//*****************************************************************************
ND::TTPCMMHoriGapMerge::TTPCMMHoriGapMerge(void) : ND::TTPCSeed2CluMatchMerge(){
  fRun = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.RunMMHorizontalGapMerging");
  if (!fRun)
    std::cout<<"TRexRecon WARNING: MM horizontal gap matching and merging disabled"<<std::endl;

  fMinDistForMatchP2J = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.MMHoriGapMerge.MinDistForMatchPathToJunction");
  fMinDistForMatchP2P = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.MMHoriGapMerge.MinDistForMatchPathToPath");
  fNbWantedRows = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.MMHoriGapMerge.NbRowsToGap");

  fAlgo = kMMHoriGapMerge;
}



//*****************************************************************************
void ND::TTPCMMHoriGapMerge::Process(ND::TReconObjectContainer *inputPatterns, ND::TReconObjectContainer *mergedPatterns){
  if ( !fRun){
    for (ND::TReconObjectContainer::iterator pattit = inputPatterns->begin(); pattit != inputPatterns->end(); pattit++) {
      mergedPatterns->push_back(*pattit);
    }
    return;
  }

  if ( ND::tpcDebug().MMHoriGapMerge(DB_INFO))
    std::cout<<"TTPCMMHoriGapMerge::Process"<<std::endl;

  // To speed things up, list patterns that start or end in each MM, of each RP, of each TPC.
  std::vector< ND::THandle<ND::TReconBase> > rawPatterns[3][2][12];
  for (ND::TReconObjectContainer::iterator pattit = inputPatterns->begin(); pattit != inputPatterns->end(); pattit++) {
    ND::THandle<ND::TTPCPattern> pattern = *pattit;
    int tpc = -1;
    int rp = -1;
    std::map<int,int> mmmap;
    FindLocation(pattern, tpc, rp, mmmap);

    pattern->SetUsable(true);
    for (std::map<int,int>::iterator MMit = mmmap.begin(); MMit != mmmap.end(); MMit++){
      rawPatterns[tpc][rp][MMit->first].push_back(pattern);
    }

    if ( ND::tpcDebug().MMHoriGapMerge(DB_VERBOSE)){
      std::cout<<" Input pattern: TPC"<<(tpc+1)<<", RP"<<rp<<", is in MM = ";
      for (std::map<int,int>::iterator MMit = mmmap.begin(); MMit != mmmap.end(); MMit++)
        std::cout<<(1+MMit->first)<<"  ";
      std::cout<<std::endl;
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

  // Loop through the TPCs and Readout planes.
  //   Loop through the patterns and find the candidates for matching and propagating.


// Define zones:
//    Here use MM = [1;12]
//    Top of MM    => Zone = MM
//    Bottom of MM => Zone = -MM
//    ZoneB = -1 * (ZoneA + (-1 * sign(ZoneA) + MM)

  std::vector< ND::THandle<ND::TReconBase> > resultPatterns;
  fMaxNbMatchCand = 0;
  fMaxNbPropagCand = 0;

  for( int tpc = 0; tpc < 3; tpc++){
    for( int rp = 0; rp < 2; rp++){
      for( int mm = 0; mm < 12; mm++){
        // Always good to short cut if there is nothing to propagate
        if ( rawPatterns[tpc][rp][mm].size() == 0)
          continue;

        for( int bt = -1; bt < 2; bt+=2){
          // ZoneA =  2 => ZoneB = -1
          // ZoneA = -2 => ZoneB =  3
          int ZoneA =    bt*(mm+1);
          int ZoneB = -1*bt*(mm+1-bt);

          if ( ND::tpcDebug().MMHoriGapMerge(DB_VERBOSE))
            std::cout<<" ==> Match-Merge for TPC "<<(tpc+1)<<", RP "<<rp<<", MM "<<(1+mm)<<", BT "<<bt<<", ZoneA "<<ZoneA<<", ZoneB "<<ZoneB<<std::endl;

          // Skip zones at top and bottom of the TPC
          if ( ZoneA == 1 || ZoneA == -6 || ZoneA == 7 || ZoneA == -12 ){
            // Save the patterns that are in 
            for (ND::TReconObjectContainer::iterator pattit = rawPatterns[tpc][rp][mm].begin(); pattit != rawPatterns[tpc][rp][mm].end(); pattit++) {
              resultPatterns.push_back(*pattit);
            }
            continue;
          }

          // ---------- Propagate the ZoneB paths to ZoneA patterns
          // Fill the containers for the downstream paths that I will propagate to upstream patterns
          for (ND::TReconObjectContainer::iterator pattit = rawPatterns[tpc][rp][mm-bt].begin(); pattit != rawPatterns[tpc][rp][mm-bt].end(); pattit++) {
            ND::THandle<ND::TTPCPattern> pattern = *pattit;
            FindPropagCandidates(pattern, ZoneB);
          }
          fMaxNbPropagCand = std::max(fNbPropagCand, fMaxNbPropagCand);
          // Process the ZoneA patterns
          for (ND::TReconObjectContainer::iterator pattit = rawPatterns[tpc][rp][mm].begin(); pattit != rawPatterns[tpc][rp][mm].end(); pattit++) {
            ND::THandle<ND::TTPCPattern> origPattern = *pattit;
            ND::THandle<ND::TTPCPattern> pattern = *pattit;
            if ( FindMatchCandidates(pattern, ZoneA)){

              fMaxNbMatchCand = std::max(fNbMatchCand, fMaxNbMatchCand);
              for (int mtcd = 0; mtcd < fNbMatchCand; mtcd++){
                FindMatches(mtcd, pattern);
                // A match was found
                if (fMatchCand[mtcd].MatchedPathId > -1){
                  // if there is a match, extend (or recreate) the upstream pattern. Mark the downstream as "used"
                  ND::THandle<ND::TTPCPattern> newPattern;
                  newPattern = MergePatterns(pattern, mtcd);

                  // If the propagCandidate is a curving back track, this new pattern may
                  // still have another state that can be propagated.
                  if ( newPattern){
                    FindPropagCandidates(newPattern, ZoneB);
                    CleanUpMatchCand();
                    FindMatchCandidates(newPattern, ZoneA);
                    pattern = newPattern;
                  }
                  fMaxNbPropagCand = std::max(fNbPropagCand, fMaxNbPropagCand);
                }
              }
            }
            // Some other parts of the this pattern might need
            // match and merge in other MMs.
            if (origPattern != pattern){
              std::map<int,int> mmmap;
              FindLocation(pattern, tpc, rp, mmmap);
              pattern->SetUsable(true);
              for (std::map<int,int>::iterator MMit = mmmap.begin(); MMit != mmmap.end(); MMit++){
                // We are already taking care of this MM.
                if ( MMit->first == (mm+1))
                  continue;
                rawPatterns[tpc][rp][MMit->first].push_back(pattern);
              }
            }
            resultPatterns.push_back(pattern);
            CleanUpMatchCand();
          }

          // ---------- Propagate the ZoneA paths to ZoneB patterns
          // Clear the containers
          fNbMatchCand = fMaxNbMatchCand;
          fNbPropagCand = fMaxNbPropagCand;
          CleanUpPropagCand();
        }
      }
    }
  }
  for (ND::TReconObjectContainer::iterator pattit = resultPatterns.begin(); pattit != resultPatterns.end(); pattit++) {
    ND::THandle<ND::TTPCPattern> pattern = *pattit;
    if ( pattern->IsUsable()){
      mergedPatterns->push_back(pattern);
      pattern->SetUsable(false);
    }
  }

  // Clean up to avoid memory leaks and so on.
  fNbMatchCand = fMaxNbMatchCand;
  fNbPropagCand = fMaxNbPropagCand;

  CleanUp();

}


//*****************************************************************************
void ND::TTPCMMHoriGapMerge::FindLocation(ND::THandle<ND::TTPCPattern> pattern, int &tpc, int &rp, std::map<int,int> &mm){
  // Loop over constituents
  ND::TReconObjectContainer::iterator constit = pattern->GetConstituents()->begin();
  for (constit = pattern->GetConstituents()->begin(); constit != pattern->GetConstituents()->end(); constit++) {
    ND::THandle<ND::TTPCPath> path = *constit;
    if (!path)
      continue;
    // Loop over first and last cluster for a path
    ND::THandle<ND::TTPCHVCluster> HVClu[2];
    HVClu[0] = *(path->GetHits()->begin());
    HVClu[1] = *(path->GetHits()->rbegin());
    for (int fl = 0; fl < 2; fl++){
      if ( !HVClu[fl])
        continue;
      ND::THandle<ND::TTPCHitPad> hitPad = *(HVClu[fl]->GetHits().begin());
      if (!hitPad){
        // TODO: Proper exception
        throw;
      }
      int geotpc, geohalf, geomm, geopad;
      ND::TTPCGeom& tpcGeom = ND::TGeomInfo::TPC();
      tpcGeom.GetGeometryInfo(hitPad->GetGeomId(), geotpc, geohalf, geomm, geopad);
      tpc = geotpc;
      rp = geohalf;
      mm[geomm]++;
    }
  }

}





//*****************************************************************************
bool ND::TTPCMMHoriGapMerge::IsClusterUsable(ND::THandle<ND::TTPCHVCluster> Cluster){
  int prevMM = -1;
  bool cluOn2MM = false;
  for (ND::THitSelection::const_iterator tHit = Cluster->GetHits().begin(); tHit != Cluster->GetHits().end(); tHit++){
    int mm = ND::TGeomInfo::TPC().GeomIdToPad((*tHit)->GetGeomId());
    if ( prevMM < 0){
      prevMM = mm;
      continue;
    }
    cluOn2MM = (prevMM != mm);
    if ( cluOn2MM )
      break;
  }

  return !cluOn2MM;
}



//*****************************************************************************
bool ND::TTPCMMHoriGapMerge::HitSelected(ND::THandle<ND::TTPCHitPad> hitPad, int &Zone){
  int mm  = ND::TGeomInfo::TPC().GeomIdToMM(hitPad->GetGeomId());
  if ( (mm+1) != abs(Zone))
    return false;

  int pad = ND::TGeomInfo::TPC().GeomIdToPad(hitPad->GetGeomId());
  int row = ND::TGeomInfo::TPC().PadToRow(pad);

  // Positive zone means top of the MM.
  // The fun of MM geometry:
  //  - for MMs 1 to 6 the top row is 0
  //  - for MMs 7 to 12 the top row is 47

  bool Select = ( (  Zone > 0 && abs(Zone) < 7 && row < (fNbWantedRows+1) ) ||
                  (  Zone > 0 && abs(Zone) > 6 && row > (47-fNbWantedRows) ) ||
                  (  Zone < 0 && abs(Zone) < 7 && row > (47-fNbWantedRows) ) ||
                  (  Zone < 0 && abs(Zone) > 6 && row < (fNbWantedRows+1) ) );

  return Select;
}
