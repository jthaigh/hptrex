#include "TTPCCathCrosserMerge.hxx"
#include "TTPCJunction.hxx"
#include "TTPCCalibration.hxx"
#include "TTPCDebug.hxx"
#include "TTPCUtils.hxx"

#include <THandle.hxx>
#include <TGeomInfo.hxx>
#include <TOARuntimeParameters.hxx>
#include <TIntegerDatum.hxx>



//*****************************************************************************
ND::TTPCCathCrosserMerge::TTPCCathCrosserMerge(ND::TTPCT0Finder *T0Finder){
  fRun = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.RunCathCrosserMerging");
  if (!fRun)
    std::cout<<"TRexRecon WARNING: cathode crossers matching and merging disabled"<<std::endl;

  fT0Finder = T0Finder;

  fCathSurf[0][0] = ND::gman().GetSetup().surface("/t2k_1/OA_0/Magnet_0/Basket_0/Tracker_0/TPC1_0/GasGap_0/Drift_0/CentralCathode_0/S2");
  fCathSurf[0][1] = ND::gman().GetSetup().surface("/t2k_1/OA_0/Magnet_0/Basket_0/Tracker_0/TPC1_0/GasGap_0/Drift_0/CentralCathode_0/S1");
  fCathSurf[1][0] = ND::gman().GetSetup().surface("/t2k_1/OA_0/Magnet_0/Basket_0/Tracker_0/TPC2_0/GasGap_0/Drift_0/CentralCathode_0/S2");
  fCathSurf[1][1] = ND::gman().GetSetup().surface("/t2k_1/OA_0/Magnet_0/Basket_0/Tracker_0/TPC2_0/GasGap_0/Drift_0/CentralCathode_0/S1");
  fCathSurf[2][0] = ND::gman().GetSetup().surface("/t2k_1/OA_0/Magnet_0/Basket_0/Tracker_0/TPC3_0/GasGap_0/Drift_0/CentralCathode_0/S2");
  fCathSurf[2][1] = ND::gman().GetSetup().surface("/t2k_1/OA_0/Magnet_0/Basket_0/Tracker_0/TPC3_0/GasGap_0/Drift_0/CentralCathode_0/S1");

  fCathSurf[2][0] = ND::gman().GetSetup().surface("/t2k_1/OA_0/Magnet_0/Basket_0/Tracker_0/TPC3_0/GasGap_0/Drift_0/CentralCathode_0/S1");
  fCathSurf[2][1] = ND::gman().GetSetup().surface("/t2k_1/OA_0/Magnet_0/Basket_0/Tracker_0/TPC3_0/GasGap_0/Drift_0/CentralCathode_0/S2");
////
////  fRangeInZ[0][0] = ND::gman().GetSetup().surface("/t2k_1/OA_0/Magnet_0/Basket_0/Tracker_0/TPC1_0/GasGap_0/Drift_0/Half_0/MM_0/S5").position()[2];
////  fRangeInZ[0][1] = ND::gman().GetSetup().surface("/t2k_1/OA_0/Magnet_0/Basket_0/Tracker_0/TPC1_0/GasGap_0/Drift_0/Half_0/MM_6/S5").position()[2];
////  fRangeInZ[1][0] = ND::gman().GetSetup().surface("/t2k_1/OA_0/Magnet_0/Basket_0/Tracker_0/TPC2_0/GasGap_0/Drift_0/Half_0/MM_0/S5").position()[2];
////  fRangeInZ[1][1] = ND::gman().GetSetup().surface("/t2k_1/OA_0/Magnet_0/Basket_0/Tracker_0/TPC2_0/GasGap_0/Drift_0/Half_0/MM_6/S5").position()[2];
////  fRangeInZ[2][0] = ND::gman().GetSetup().surface("/t2k_1/OA_0/Magnet_0/Basket_0/Tracker_0/TPC3_0/GasGap_0/Drift_0/Half_0/MM_0/S5").position()[2];
////  fRangeInZ[2][1] = ND::gman().GetSetup().surface("/t2k_1/OA_0/Magnet_0/Basket_0/Tracker_0/TPC3_0/GasGap_0/Drift_0/Half_0/MM_6/S5").position()[2];
////
////  // This is not for a high preceision test. One TPC is good enough for all tests.
  fRangeInY[0] = ND::gman().GetSetup().surface("/t2k_1/OA_0/Magnet_0/Basket_0/Tracker_0/TPC2_0/GasGap_0/Drift_0/S3").position()[1];
  fRangeInY[1] = ND::gman().GetSetup().surface("/t2k_1/OA_0/Magnet_0/Basket_0/Tracker_0/TPC2_0/GasGap_0/Drift_0/S4").position()[1];

  fMaxReduChi2Match = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.CathodeCrosserMerge.MaxReduChi2Match");
  fMaxFirstResidual = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.CathodeCrosserMerge.MaxFirstResidual");

  fNbMatchCand = MAXNBCATHMATCHCANDIDATE;
  fNbPropagCand = MAXNBCATHPROPAGCANDIDATE;
  CleanUp();
}



//*****************************************************************************
void ND::TTPCCathCrosserMerge::Process(ND::TReconObjectContainer *inputPatterns, ND::TReconObjectContainer *mergedPatterns){
  if ( !fRun){
    for (ND::TReconObjectContainer::iterator pattit = inputPatterns->begin(); pattit != inputPatterns->end(); pattit++) {
      mergedPatterns->push_back(*pattit);
    }
    return;
  }

    

  if ( ND::tpcDebug().CathCrosserMerge(DB_INFO))
    std::cout<<"TTPCCathCrosserMerge::Process"<<std::endl;

  std::vector< ND::THandle<ND::TReconBase> > rawPatterns[3][2];
  for (ND::TReconObjectContainer::iterator pattit = inputPatterns->begin(); pattit != inputPatterns->end(); pattit++) {
    ND::THandle<ND::TTPCPattern> pattern = *pattit;
    int tpc = -1;
    int rp = -1;
    FindLocation(pattern, tpc, rp);
    pattern->SetUsable(true);
    rawPatterns[tpc][rp].push_back(pattern);
    if ( ND::tpcDebug().CathCrosserMerge(DB_VERBOSE)){
      std::cout<<" Input pattern: TPC"<<(tpc+1)<<", RP"<<rp<<std::endl;
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

  TTPCRecPackUtils::SavedParam RPParam = TTPCRecPackUtils::InitForPropagationInTPC();
  // Enable energy loss correction
  ND::rpman("TREx").model_svc().enable_correction(RP::particle_helix, RP::eloss, true);
  ND::rpman("TREx").model_svc().enable_noiser(RP::particle_helix, RP::eloss, true);
  ND::rpman("TREx").model_svc().enable_noiser(RP::particle_helix, RP::ms, true);
  ND::rpman("TREx").model_svc().model().intersector().set_length_sign(1);

  fNbMatchCand = 0;
  fNbPropagCand = 0;
  fMaxNbMatchCand = 0;
  fMaxNbPropagCand = 0;
  // Loop through the TPCs and Readout planes.
  //   Loop through the patterns and find the candidates for matching and propagating.

  for( int tpc = 0; tpc < 3; tpc++){
    if ( ND::tpcDebug().CathCrosserMerge(DB_VERBOSE))
      std::cout<<" ==> Match patterns in TPC "<<(tpc+1)<<std::endl;
    std::vector< ND::THandle<ND::TReconBase> > resultPatterns;
    // ---------- Propagate from RP1
    int rpMatch = 0;
    int rpPropa = 1;
    // Clear the containers
    // CleanUp();
    // Fill the containers for the RP1 paths that I will propagate to RP0 patterns
    for (ND::TReconObjectContainer::iterator pattit = rawPatterns[tpc][rpPropa].begin(); pattit != rawPatterns[tpc][rpPropa].end(); pattit++) {
      ND::THandle<ND::TTPCPattern> pattern = *pattit;
      FindPropagCandidates(pattern, tpc, rpMatch);
    }

    // Process the RP0 patterns
    if ( ND::tpcDebug().CathCrosserMerge(DB_VERBOSE))
      std::cout<<"   Match to the RP0 patterns"<<std::endl;
    MatchPatternsAWithB(rawPatterns[tpc][rpMatch], rpMatch, rpPropa, resultPatterns);

    // ---------- Propagate from RP0
    rpMatch = 1;
    rpPropa = 0;
    // Clear the containers
    fNbMatchCand = fMaxNbMatchCand;
    fNbPropagCand = fMaxNbPropagCand;
    CleanUp();
    // Fill the containers for the RP0 paths (from resultPatterns) that I will propagate to RP1 patterns
    for (ND::TReconObjectContainer::iterator pattit = resultPatterns.begin(); pattit != resultPatterns.end(); pattit++) {
      ND::THandle<ND::TTPCPattern> pattern = *pattit;
      pattern->SetUsable(true);
      FindPropagCandidates(pattern, tpc, rpMatch);
    }

    // Process the RP1 patterns
    if ( ND::tpcDebug().CathCrosserMerge(DB_VERBOSE))
      std::cout<<"   Match to the RP1 patterns"<<std::endl;
    MatchPatternsAWithB(rawPatterns[tpc][rpMatch], rpMatch, rpPropa, resultPatterns);

    for (ND::TReconObjectContainer::iterator pattit = resultPatterns.begin(); pattit != resultPatterns.end(); pattit++) {
      ND::THandle<ND::TTPCPattern> pattern = *pattit;
      if ( ! pattern->IsUsable())
        continue;
      mergedPatterns->push_back(pattern);
    }
    // Clear the containers
    fNbMatchCand = fMaxNbMatchCand;
    fNbPropagCand = fMaxNbPropagCand;
    CleanUp();
    
  }
  
  // Clean up to avoid memory leaks and so on.
  fNbMatchCand = fMaxNbMatchCand;
  fNbPropagCand = fMaxNbPropagCand;
  CleanUp();
  TTPCRecPackUtils::ResetAfterPropagationInTPC(RPParam);
}


//*****************************************************************************
void ND::TTPCCathCrosserMerge::FindLocation(ND::THandle<ND::TTPCPattern> pattern, int &tpc, int &rp){
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
}


//*****************************************************************************
void ND::TTPCCathCrosserMerge::FindPropagCandidates(ND::THandle<ND::TTPCPattern> pattern, int tpc, int matchRP){
  if ( ! pattern->IsUsable())
    return;
  if ( !(fNbPropagCand < MAXNBCATHPROPAGCANDIDATE))
    return;

  if ( ND::tpcDebug().CathCrosserMerge(DB_VERBOSE))
    std::cout<<" ---- FindPropagCandidates"<<std::endl;
  for (ND::TReconObjectContainer::iterator constit = pattern->GetConstituents()->begin(); constit != pattern->GetConstituents()->end(); constit++) {
    ND::THandle<ND::TTPCPath> path = (*constit);

    // Here prevent creation of too many propag candidates.
    if ( !(fNbPropagCand < MAXNBCATHPROPAGCANDIDATE))
      break;
    if (! path)
      continue;
    if (!path->HasFitState())
      continue;

////    ND::tman("TREx").SetVerbosity("07000");
////    ND::rpman("TREx").navigation_svc().navigator(RP::particle_helix).set_unique_surface(false);
    State frontState = path->GetFrontFitState();
    State backState  = path->GetBackFitState();
    if(frontState.name(RP::representation) != RP::pos_dir_curv)
      RP::rep().convert(frontState, RP::pos_dir_curv);
    // Check where the path state is pointing to.
    // If it is not pointing toward the cathode, revert it
////    std::cout<<" >=>>=>> "<<matchRP<<"   "<<frontState.vector()[3]<<std::endl;
    if ( (frontState.vector()[3] > 0 && matchRP == 0) || (frontState.vector()[3] < 0 && matchRP == 1)){
      ND::tman().ReverseStateSenseAndCharge(frontState);
      ND::tman().ReverseStateSenseAndCharge(backState);
    }

    // Propagate to the other side of the cathode and verify if the propagation point is within YZ boundaries
    // of the drift volume.
    // Surface surf = ND::rpman().geometry_svc().setup().
////    double Length;
////    std::cout<<" ==>>=>> "<<fCathSurf[tpc][matchRP].position()[0]<<"   "<<propagState.vector()[0]<<"   "<<propagState.vector()[1]<<"   "<<propagState.vector()[2]<<"   "<<propagState.vector()[3]<<std::endl;
////    bool ok = ND::rpman("TREx").navigation_svc().propagate(fCathSurf[tpc][matchRP], propagState, Length);                                                       
////    std::cout<<" ==>>=>> "<<fCathSurf[tpc][matchRP].position()[0]<<"   "<<propagState.vector()[0]<<"   "<<propagState.vector()[1]<<"   "<<propagState.vector()[2]<<"   "<<propagState.vector()[3]<<"   "<<ok<<std::endl;
////    if (!ok)
////      continue;
////
////    // Is the propagated state within the TPC boundaries ?
////    std::cout<<" ===>>>> "<<propagState.vector()[1]<<" < "<<fRangeInY[0]<<" || "<<propagState.vector()[1]<<" > "<<fRangeInY[1]<<std::endl;
////    if (propagState.vector()[1] < fRangeInY[0] || propagState.vector()[1] > fRangeInY[1])
////      continue;
    // TODO: Probably remove this since I don't necessarily have the T0
    //// if (propagState.vector()[2] < fRangeInZ[tpc][0] || propagState.vector()[2] > fRangeInZ[tpc][1])
    ////   continue;

    // If it is, store the state in PropagCandidate
    // Convert X position into time too because we want to match in the TYZ space, not XYZ space in case one of the segments doesn't have a T0, or not the right one.
    //// propagState.vector()[0] = (propagState.vector()[0] / ND::tpcCalibration().GetDriftVelocity()) + path->GetT0() + ND::tpcCalibration().GetTimeOffset();

    fPropagCand[fNbPropagCand].Pattern     = pattern;
    fPropagCand[fNbPropagCand].Path        = path;
    fPropagCand[fNbPropagCand].propagState[0] = frontState;
    fPropagCand[fNbPropagCand].propagState[1] = backState;

    fNbPropagCand++;
    if ( !(fNbPropagCand < MAXNBCATHPROPAGCANDIDATE))
      return;
  }
//// ND::tman("TREx").SetVerbosity("00000");
}


//*****************************************************************************
void ND::TTPCCathCrosserMerge::CleanUp() {
  CleanUpMatchCand();
  CleanUpPropagCand();
}


//*****************************************************************************
void ND::TTPCCathCrosserMerge::CleanUpMatchCand() {
  for (int mtcd = 0; mtcd < fNbMatchCand; mtcd++){
    fMatchCand[mtcd].Path = ND::THandle<ND::TTPCPath>();
    fMatchCand[mtcd].ClustersToMatch = ND::THandle<ND::THitSelection>();
    fMatchCand[mtcd].MatchedPathId = -1;
    fMatchCand[mtcd].ReduChi2 = fMaxReduChi2Match;
  }
  fNbMatchCand = 0;
}


//*****************************************************************************
void ND::TTPCCathCrosserMerge::CleanUpPropagCand() {
  for (int ppcd = 0; ppcd < fNbPropagCand; ppcd++){
    fPropagCand[ppcd].Pattern = ND::THandle<ND::TTPCPattern>();
    fPropagCand[ppcd].Path = ND::THandle<ND::TTPCPath>();
    for (int i = 0; i < 2; i++){
      HyperVectorObject &tmpHVO = fPropagCand[ppcd].propagState[i];
      NamedObject &tmpNO = fPropagCand[ppcd].propagState[i];
      tmpHVO.clear();
      tmpNO.clear();
    }
  }
  fNbPropagCand = 0;
  
}




//*****************************************************************************
void ND::TTPCCathCrosserMerge::MatchPatternsAWithB(std::vector< ND::THandle<ND::TReconBase> > &rawPat, int zoneA, int zoneB, std::vector< ND::THandle<ND::TReconBase> > &resPat){
  std::vector< ND::THandle<ND::TReconBase> > tmpPat;

  // FindPropagCandidates should have been called just before this method.
  // Best time to check the max number of propag candidates.
  fMaxNbPropagCand = std::max(fNbPropagCand, fMaxNbPropagCand);
  for (ND::TReconObjectContainer::iterator pattit = rawPat.begin(); pattit != rawPat.end(); pattit++) {
    ND::THandle<ND::TTPCPattern> pattern = *pattit;
    if ( FindMatchCandidates(pattern, zoneA)){
      fMaxNbMatchCand = std::max(fNbMatchCand, fMaxNbMatchCand);
      for (int mtcd = 0; mtcd < fNbMatchCand; mtcd++){
        FindMatches(mtcd, pattern);
        // A match was found
        if (fMatchCand[mtcd].MatchedPathId > -1){
          // if there is a match, extend (or recreate) the matched pattern. Mark the propagated path as "used"
          ND::THandle<ND::TTPCPattern> newPattern;
          newPattern = MergePatterns(pattern, mtcd);
          // A pattern can have multiple paths crossing the cathode
          if ( newPattern){
            CleanUpMatchCand();
            FindMatchCandidates(newPattern, zoneA);
            mtcd = -1;
            pattern = newPattern;
          }
        }
      }
      CleanUpMatchCand();
    }
    tmpPat.push_back(pattern);
  }
  for (ND::TReconObjectContainer::iterator pattit = tmpPat.begin(); pattit != tmpPat.end(); pattit++) {
    ND::THandle<ND::TTPCPattern> pattern = *pattit;
    if ( ! pattern->IsUsable()){
      continue;
    }
    resPat.push_back(pattern);
  }
  tmpPat.clear();
}



//*****************************************************************************
bool ND::TTPCCathCrosserMerge::FindMatchCandidates(ND::THandle<ND::TTPCPattern> pattern, int zone){
  if ( ND::tpcDebug().CathCrosserMerge(DB_VERBOSE))
    std::cout<<" ---- FindMatchCandidates"<<std::endl;
  if ( ! pattern->IsUsable())
    return false;

  // For each pattern find the paths and junctions at the edge of the MM. Fill simple containers without pattern in it.
  bool FoundCandidate = false;
  for (ND::TReconObjectContainer::iterator constit = pattern->GetConstituents()->begin(); constit != pattern->GetConstituents()->end(); constit++) {
    // Here prevent creation of too many match candidates.
    if ( !(fNbMatchCand < MAXNBCATHMATCHCANDIDATE))
      break;
    ND::THandle<ND::TTPCPath> path = (*constit);
    if (!path)
      continue;
    if (path->GetIsXPath())
      continue;
    // Skip the cathode crossers already merged
    if (path->IsCathodeCrosser())
      continue;

    // TODO: Make a selection based on the location of the junction
    // with respect to the cathode.
    bool selected = true;

    if (selected){
      FoundCandidate = true;
      // These residuals define the minimum distance between propagated seed and
      // a cluster needed to consider that we have a match
      fMatchCand[fNbMatchCand].Path = path;
      fMatchCand[fNbMatchCand].ReduChi2 = fMaxReduChi2Match;

      ND::THandle<ND::TTPCHVCluster> FirstClu = *(path->GetHits()->begin());
      ND::THandle<ND::TTPCHVCluster> LastClu  = *(path->GetHits()->rbegin());
      bool reorder = false;
      if (zone == 0){
        reorder = FirstClu->CalibX() < LastClu->CalibX();
      } else {
        reorder = FirstClu->CalibX() > LastClu->CalibX();
      }

      ND::THandle<ND::THitSelection> Clusters;
      if( reorder ){
        Clusters = ND::THandle<ND::THitSelection>(new ND::THitSelection());
        for (ND::THitSelection::const_reverse_iterator Clu = path->GetHits()->rbegin(); Clu != path->GetHits()->rend(); Clu++) {
          Clusters->push_back(*Clu);
        }
      } else {
        Clusters = path->GetHits();
      }
      fMatchCand[fNbMatchCand].ClustersToMatch = Clusters;

      fNbMatchCand++;

      
    }

  }

  return FoundCandidate;
}


//*****************************************************************************
void ND::TTPCCathCrosserMerge::FindMatches(int matchId, ND::THandle<ND::TTPCPattern> MatchCandPat){
  // Loop over the MatchCandidates
  //   Loop over the fPathToPropag
  //     Propagate fPathToPropag to the cluster or pads
  //      if distance less than cutoff (10cm ?) then compare with previous residuals
  //      if residuals smaller than saved one, then save the new one.

  // We have to force the sense of the propagation here for
  // for curving back tracks in particular.

  // If hit to match is HVCluster => easy
  if ( ND::tpcDebug().CathCrosserMerge(DB_VERBOSE))
    std::cout<<" ---- FindMatches to HV cluster from a path"<<std::endl;

  if ( ND::tpcDebug().CathCrosserMerge(DB_VERBOSE)){
    ND::THandle<ND::TTPCPath> tmpPath = fMatchCand[matchId].Path;
    if (tmpPath)
      std::cout<<"    -> MatchCand, path Id "<<tmpPath->GetId()<<std::endl;
  }

  for (int ppcd = 0; ppcd < fNbPropagCand; ppcd++){
    if ( ! fPropagCand[ppcd].Pattern->IsUsable())
      continue;
    // Don't match a pattern with itself.
    // Important when the pattern is on top of a MM horizontal gap.
    if ( MatchCandPat == fPropagCand[ppcd].Pattern)
      continue;
    if ( ND::tpcDebug().CathCrosserMerge(DB_VERBOSE)){
      std::cout<<"    -> PropagCand, path Id "<<fPropagCand[ppcd].Path->GetId()<<std::endl;
    }

    ND::rpman("TREx").model_svc().model().intersector().set_length_sign(1);
    State propagState;
    bool StateFound = false;
    double Distance = 9999.9;
    for (int i = 0; i < 2; i++){
      propagState = fPropagCand[ppcd].propagState[i];
      // std::cout<<" ===>>> State     "<<propagState.vector()[0]<<"  "<<propagState.vector()[1]<<"  "<<propagState.vector()[2]<<std::endl;
      ND::THandle<ND::TTPCHVCluster> TargetCluster = *(fMatchCand[matchId].ClustersToMatch->begin());
      TVector3 clusterPosition = TargetCluster->GetCalibratedPosition();

      if (!TTPCRecPackUtils::PropagateToHVCluster(TargetCluster, propagState)){
        if ( ND::tpcDebug().CathCrosserMerge(DB_VVERBOSE))
          std::cout<<"       RecPack propagation to first cluster FAILED for propagation state "<<(i+1)<<std::endl;
        continue;
      }
      // DO NOT MATCH X because we potentially don't have the proper or any T0 for one or both segments.
      TVector3 propagPosition(clusterPosition.X(), propagState.vector()[1], propagState.vector()[2]);
      // std::cout<<" ===>>> Cluster   "<<clusterPosition.X()<<"  "<<clusterPosition.Y()<<"  "<<clusterPosition.Z()<<std::endl;
      // std::cout<<" ===>>> Propag    "<<propagPosition.X()<<"  "<<propagPosition.Y()<<"  "<<propagPosition.Z()<<std::endl;
      Distance = (propagPosition - clusterPosition).Mag();
      if (Distance < fMaxFirstResidual){
        StateFound = true;
        break;
      }
    }
    if (StateFound){
      if ( ND::tpcDebug().CathCrosserMerge(DB_VERBOSE))
        std::cout<<"       First cluster check SUCCESSFUL with distance = "<<Distance<<std::endl;
    } else {
      if ( ND::tpcDebug().CathCrosserMerge(DB_VERBOSE))
        std::cout<<"       First cluster check FAILED with distance = "<<Distance<<std::endl;
      continue;
    }

    // Time for the chi2 check !
    double Chi2_dir[3];
    int NDOF_dir[3]; 
   
    for (int i = 0; i < 3; i++){
      Chi2_dir[i] = 0.0;
      NDOF_dir[i] = 0;
    }
    double Chi2 = 0.0;
    int NDOF = 0;
    
    //// ND::tman("TREx").SetVerbosity("05000");
    // Use the First ClustersToMatch to know if the clusters should be reordered
    ND::THandle<ND::THitSelection> Clusters = fMatchCand[matchId].ClustersToMatch;

    State prevState = propagState;
    for (ND::THitSelection::iterator tmpClu = Clusters->begin(); tmpClu != Clusters->end(); tmpClu++) {
      ND::THandle<ND::TTPCHVCluster> Clu = *tmpClu;

      if ( ND::tpcDebug().CathCrosserMerge(DB_VVERBOSE)){
        std::cout<<"    -> Propagate from "<<propagState.vector()[0]<<"  "<<propagState.vector()[1]<<"  "<<propagState.vector()[2]<<std::endl;
        std::cout<<"       to cluster at  "<<Clu->GetPosition().X()<<"  "<<Clu->GetPosition().Y()<<"  "<<Clu->GetPosition().Z()<<std::endl;
      }
      double length;
      if (!TTPCRecPackUtils::FullPropagateToHVCluster(Clu, propagState, length)){
        if ( ND::tpcDebug().CathCrosserMerge(DB_VVERBOSE))
          std::cout<<"       Propagation FAILED"<<std::endl;
        continue;
      }
      if ( length > 2000.){
        if ( ND::tpcDebug().CathCrosserMerge(DB_VVERBOSE))
          std::cout<<"       Propagation length = "<<length<<"mm => TOO LONG ! Skip propagation."<<std::endl;
        propagState = prevState; // Go back to the previous iteration
        continue;
      }
      if ( ND::tpcDebug().CathCrosserMerge(DB_VVERBOSE))
        std::cout<<"    -> Propagated to  "<<propagState.vector()[0]<<"  "<<propagState.vector()[1]<<"  "<<propagState.vector()[2]<<std::endl;
      TTPCUtils::State2CluChi2(propagState, Clu, Chi2_dir, NDOF_dir);

      // State2CluChi2 already counts according to the cluster orientation
      Chi2 += Chi2_dir[1] + Chi2_dir[2];
      NDOF += NDOF_dir[1] + NDOF_dir[2];

      if ( ND::tpcDebug().CathCrosserMerge(DB_VVERBOSE))
        std::cout<<"       Propagation successful with reduced chi2 = "<<(Chi2_dir[1] + Chi2_dir[2])<<std::endl;
      prevState = propagState;
    }

    NDOF -= 3; // Take into account the number of fitting parameters

    // Make sure that enough clusters were actually matched.
    // Matching a segment <10 clusters does matter anyway.
    if (NDOF > 10) {
      if ( (Chi2/double(NDOF)) < fMatchCand[matchId].ReduChi2 ){
        fMatchCand[matchId].MatchedPathId = ppcd;
        fMatchCand[matchId].ReduChi2 = Chi2/double(NDOF);

        if ( ND::tpcDebug().CathCrosserMerge(DB_VERBOSE))
          std::cout<<"       New best match with reduced chi2 = "<<Chi2<<" / "<<NDOF<<" = "<<(Chi2/double(NDOF))<<"   ppcd = "<<ppcd<<std::endl;
      } else {
        if ( ND::tpcDebug().CathCrosserMerge(DB_VERBOSE))
          std::cout<<"       Match FAILED with reduced chi2 = "<<Chi2<<" / "<<NDOF<<" = "<<(Chi2/double(NDOF))<<"   ppcd = "<<ppcd<<std::endl;
      }
    } else {
      if ( ND::tpcDebug().CathCrosserMerge(DB_VERBOSE))
        std::cout<<"       Match FAILED with NDOF = "<<NDOF<<"   ppcd = "<<ppcd<<std::endl;
    }
  }

//// ND::tman("TREx").SetVerbosity("00000");

}



//*****************************************************************************
ND::THandle<ND::TTPCPattern> ND::TTPCCathCrosserMerge::MergePatterns(ND::THandle<ND::TTPCPattern> origPattern, int matchId){
  if ( ND::tpcDebug().CathCrosserMerge(DB_VERBOSE))
    std::cout<<" ---- MergePatterns"<<std::endl;
  ND::THandle<ND::TTPCPattern> newPattern;

  int ppcd = fMatchCand[matchId].MatchedPathId;
  ND::THandle<ND::TTPCPath> propPath = fPropagCand[ppcd].Path;
  ND::THandle<ND::TTPCPath> matchPath = fMatchCand[matchId].Path;
  ///////////////// path + path
  if ( matchPath) {

    // ====> 1) Paths merging
    ND::THandle<ND::TTPCPath> newPath = TTPCUtils::MergePaths( matchPath, propPath);
    newPath->SetCathodeCrosser(true);
    newPath->AddConstituent(matchPath);
    newPath->AddConstituent(propPath);
    TTPCT0 newT0;
    
    // In order to recover bad T0s, it's better to just always recalculate the T0
    // from the hits at the cathode even if we already have T0 values.
    fT0Finder->FindCathodeCrosserT0(matchPath, propPath, newT0);
    origPattern->SetT0(newT0);
    fPropagCand[ppcd].Pattern->SetT0(newT0);

    // Make sure we don't reuse these patterns for another merging.
    origPattern->SetUsable(false);
    fPropagCand[ppcd].Pattern->SetUsable(false);

    // ====> 2) Pattern merging
    if ( ND::tpcDebug().CathCrosserMerge(DB_INFO))
      std::cout<<"   Merging MatchCand path "<<matchPath->GetId()<<" and PropagCand path "<<propPath->GetId()<<" into the new path "<<newPath->GetId()<<std::endl;

    //   IF two single tracks, DONE HERE !
    if ( matchPath->GetNJunctionIds() == 0 && propPath->GetNJunctionIds() == 0){
      if ( ND::tpcDebug().CathCrosserMerge(DB_VERBOSE))
        std::cout<<"   -> Simple track-track merging"<<std::endl;
      newPattern = ND::THandle<ND::TTPCPattern>( new ND::TTPCPattern(newPath) );
      newPattern->SetId(ND::tpcCalibration().GetPatternId());
      newPattern->InitialSetup();

    }
    // Create a pattern from scratch with the new junction(s) connected to the merged path and the untouched junctions.
    else {
      newPattern = ND::THandle<ND::TTPCPattern>( new ND::TTPCPattern() );
      ND::THandle<ND::TTPCPattern> oldPattern[2];
      ND::THandle<ND::TTPCPath> oldPath[2];
      int nComposite = 0; // Number of input patterns with at least one junction, a.k.a. "composite".

      // Which one has junctions ?
      // The pattern with only one path can just be forgotten because its path is in the new merged path.
      if ( matchPath->GetNJunctionIds() > 0){
        oldPattern[nComposite] = origPattern;
        oldPath[nComposite] = matchPath;
        nComposite++;
      }
      if ( propPath->GetNJunctionIds() > 0){
        oldPattern[nComposite] = fPropagCand[ppcd].Pattern;
        oldPath[nComposite] = propPath;
        nComposite++;
      }

      if ( ND::tpcDebug().CathCrosserMerge(DB_VERBOSE)){
        if ( nComposite == 1 ) {
          std::cout<<"   -> track-pattern merging"<<std::endl;
        } else if ( nComposite == 2 ) {
          std::cout<<"   -> pattern-pattern merging"<<std::endl;
        } else {
          std::cout<<"   MERGING WARNING: The number of composite pattern is "<<nComposite<<" => This should NOT happen !"<<std::endl;
        }
      }
        
      // Loop through the junctions
      for (int nC = 0; nC < nComposite; nC++){
        MigratePatterns(oldPattern[nC], oldPath[nC]->GetId(), newPath, newPattern);
        // Copy and store untouched paths.
      }
      newPattern->SetId(ND::tpcCalibration().GetPatternId());
      newPattern->InitialSetup();
    }

    // properly reset the match variables
    fMatchCand[matchId].ReduChi2 = fMaxReduChi2Match;
    fMatchCand[matchId].MatchedPathId = -1;

    newPattern->SetT0(newT0);
    return newPattern;
  }

  std::cerr<<" This is a very bad day for you. The matched object is not a TTPCPath. Let's crash !"<<std::endl;
  throw;
}



//*****************************************************************************
void ND::TTPCCathCrosserMerge::MigratePatterns(ND::THandle<ND::TTPCPattern> PatternA, unsigned int PathIdA, ND::THandle<ND::TTPCPath> PathAB, ND::THandle<ND::TTPCPattern> PatternB){

  if ( ND::tpcDebug().CathCrosserMerge(DB_VVERBOSE))
    std::cout<<"   Migrate pattern "<<PatternA->GetId()<<" into the new pattern "<<PatternB->GetId()<<std::endl;
  // 1) Recreate/copy the unmerged paths
  // Use map to keep match between old and new paths
  std::map< ND::THandle<ND::TTPCPath>,ND::THandle<ND::TTPCPath> > PathCopy;
  for (ND::TReconObjectContainer::iterator constit = PatternA->GetConstituents()->begin(); constit != PatternA->GetConstituents()->end(); constit++) {
    ND::THandle<ND::TTPCPath> tmpPath = (*constit);
    if (!tmpPath)
      continue;

    // Only copy the unmerged ones
    if ( tmpPath->GetId() == PathIdA)
      continue;

    ND::THandle<ND::TTPCPath> newPath( new TTPCPath(*tmpPath));
    // If this path has not already been copied or merged,
    // we need to store the original (i.e. gas output) path as constituent.
    if (tmpPath->GetId() < TREXSTDOUTIDOFFSET){
      newPath->AddConstituent(tmpPath);
    }
    newPath->SetId(ND::tpcCalibration().GetPathId());
    newPath->ClearJunctionIds();
    PathCopy[tmpPath] = newPath;
    if ( ND::tpcDebug().CathCrosserMerge(DB_VVERBOSE))
      std::cout<<"   + Old path Id "<<tmpPath->GetId()<<" copied into new path Id "<<newPath->GetId()<<std::endl;

  }

  // Recreate all the junctions, making sure to properly
  // store the merged path instead of the original segments.
  for (ND::TReconObjectContainer::iterator constit = PatternA->GetConstituents()->begin(); constit != PatternA->GetConstituents()->end(); constit++) {
    // Find the junction that is connected to the merged path in the old pattern.
    ND::THandle<ND::TTPCJunction> JunctionA = (*constit);
    if (!JunctionA)
      continue;

    // Copy JunctionA position
    ND::THandle<ND::TTPCJunction> newJunction( new ND::TTPCJunction(JunctionA->GetPosition()) );
    // Copy JunctionA hits
    ND::THandle<ND::THitSelection> oldHits = JunctionA->GetHits();
    ND::THitSelection *newHits = new ND::THitSelection();
    for (ND::THitSelection::const_iterator tmpHit = oldHits->begin(); tmpHit != oldHits->end(); tmpHit++) {
      newHits->push_back(*tmpHit);
    }
    newJunction->AddHits(newHits);
    // Copy JunctionA paths except for merged paths
    for (ND::TReconObjectContainer::iterator pathtit = JunctionA->GetConstituents()->begin(); pathtit != JunctionA->GetConstituents()->end(); pathtit++) {
      ND::THandle<ND::TTPCPath> tmpPath = (*pathtit);
      if ( tmpPath->GetId() == PathIdA){
        newJunction->AddConstituent(PathAB);
      } else {
        newJunction->AddConstituent(PathCopy[tmpPath]);
      }
    }
    newJunction->SetId(ND::tpcCalibration().GetJunctionId());
    if ( ND::tpcDebug().CathCrosserMerge(DB_VVERBOSE))
      std::cout<<"   + Old junction Id "<<JunctionA->GetId()<<" copied into new junction Id "<<newJunction->GetId()<<std::endl;
    PatternB->AddJunction(newJunction);
  }
}


