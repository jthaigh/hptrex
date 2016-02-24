#include "TTPCSeed2CluMatchMerge.hxx"

#include "TTPCJunction.hxx"
#include "TTPCRecPackUtils.hxx"
#include "TTPCCalibration.hxx"
#include "TTPCDebug.hxx"
#include "TTPCUtils.hxx"

#include <THandle.hxx>
#include <TGeomInfo.hxx>




//*****************************************************************************
ND::TTPCSeed2CluMatchMerge::TTPCSeed2CluMatchMerge(void){

  // These default values must be overwritten in the class inheriting from this one.
  fMinDistForMatchP2P = 0.;
  fMinDistForMatchP2J = 0.;

  fTPCBits[0] = ND::TReconBase::kTPC1;
  fTPCBits[1] = ND::TReconBase::kTPC2;
  fTPCBits[2] = ND::TReconBase::kTPC3;

  fNbMatchCand  = MAXNBMATCHCANDIDATE;
  fNbPropagCand = MAXNBPROPAGCANDIDATE;

  // Loop to clear all the THandles of fMatchCand
  CleanUp();

  fNbMatchCand  = 0;
  fNbPropagCand = 0;

  // Seeding
  fNewSeeding = new ND::TTPCSeeding();

  fAlgo = kNoSeed2Clu;

  fFakeJunctionForStudies = false;
}


//*****************************************************************************
void ND::TTPCSeed2CluMatchMerge::CleanUp() {
  CleanUpMatchCand();
  CleanUpPropagCand();
}


//*****************************************************************************
void ND::TTPCSeed2CluMatchMerge::CleanUpMatchCand() {
  for (int mtcd = 0; mtcd < fNbMatchCand; mtcd++){
    fMatchCand[mtcd].PathOrJunction = ND::THandle<ND::TReconBase>();
    fMatchCand[mtcd].ClusterToMatch[0] = ND::THandle<ND::TTPCHVCluster>();
    fMatchCand[mtcd].ClusterToMatch[1] = ND::THandle<ND::TTPCHVCluster>();
    fMatchCand[mtcd].NbClusterToMatch = 0;
    fMatchCand[mtcd].JuncPosToMatch = TVector3(-999999.9,-999999.9,-999999.9);
    fMatchCand[mtcd].MatchedPathId = -1;
    fMatchCand[mtcd].MatchedClusterIdx = -1;
    fMatchCand[mtcd].Residuals = TVector3(fMinDistForMatchP2P, fMinDistForMatchP2P, fMinDistForMatchP2P);
  }
  fNbMatchCand = 0;
}


//*****************************************************************************
void ND::TTPCSeed2CluMatchMerge::CleanUpPropagCand() {
  for (int ppcd = 0; ppcd < fNbPropagCand; ppcd++){
    fPropagCand[ppcd].Pattern = ND::THandle<ND::TTPCPattern>();
    fPropagCand[ppcd].Path = ND::THandle<ND::TTPCPath>();
    HyperVectorObject & tmpHVO = fPropagCand[ppcd].propagState;
    NamedObject & tmpNO = fPropagCand[ppcd].propagState;
    tmpHVO.clear();
    tmpNO.clear();
    fPropagCand[ppcd].FrontOrBack = -1;
  }
  fNbPropagCand = 0;
  
}


//*****************************************************************************
void ND::TTPCSeed2CluMatchMerge::MatchPatternsAWithB(std::vector< ND::THandle<ND::TReconBase> > &rawPat, int zoneA, int zoneB, std::vector< ND::THandle<ND::TReconBase> > &resPat){
  std::vector< ND::THandle<ND::TReconBase> > tmpPat;

  // FindPropagCandidates should have been called just before this method.
  // Best time to check the max number of propag candidates.
  fMaxNbPropagCand = std::max(fNbPropagCand, fMaxNbPropagCand);
  for (ND::TReconObjectContainer::iterator pattit = rawPat.begin(); pattit != rawPat.end(); pattit++) {
    ND::THandle<ND::TTPCPattern> pattern = *pattit;
    if ( FindMatchCandidates(pattern, zoneA)){

      for (int mtcd = 0; mtcd < fNbMatchCand; mtcd++){
        FindMatches(mtcd, pattern);
        fMaxNbMatchCand = std::max(fNbMatchCand, fMaxNbMatchCand);
        // A match was found
        if (fMatchCand[mtcd].MatchedPathId > -1){
          // if there is a match, extend (or recreate) the upstream pattern. Mark the downstream as "used"
          ND::THandle<ND::TTPCPattern> newPattern;
          newPattern = MergePatterns(pattern, mtcd);

          // If the propagCandidate is a curving back track, this new pattern may
          // still have another state that can be propagated.
          // Also if we matched a path to a junction at the edge of the MM,
          // then maybe another path also matches to this junction so let's
          // go over this pattern all over again.
          if ( newPattern){
            FindPropagCandidates(newPattern, zoneB);
            CleanUpMatchCand();
            FindMatchCandidates(newPattern, zoneA);
            mtcd = -1;
            pattern = newPattern;
          }
          fMaxNbPropagCand = std::max(fNbPropagCand, fMaxNbPropagCand);
        }
      }
      // Clean up only the MatchCand because of curving back tracks still present as PropagCand
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
void ND::TTPCSeed2CluMatchMerge::FindPropagCandidates(ND::THandle<ND::TTPCPattern> pattern, int zone){
  if ( ! pattern->IsUsable())
    return;

  if ( ND::tpcDebug().MatchAndMerge(DB_VERBOSE))
    std::cout<<" ---- FindPropagCandidates"<<std::endl;
  for (ND::TReconObjectContainer::iterator constit = pattern->GetConstituents()->begin(); constit != pattern->GetConstituents()->end(); constit++) {
    ND::THandle<ND::TTPCPath> path = (*constit);

    // Here prevent creation of too many propag candidates.
    if ( !(fNbPropagCand < MAXNBPROPAGCANDIDATE))
      break;
    if (! path)
      continue;
    if (!path->HasSeedState())
      continue;

    int SelectedClu;
    GetPathEndsAtZoneEdge(path, SelectedClu, zone);
    if (SelectedClu > -1){
      // Need to save the first or last node state depending on the cluster
      // selected here.
      if (SelectedClu == 0 || SelectedClu == 2){
        fPropagCand[fNbPropagCand].Pattern = pattern;
        fPropagCand[fNbPropagCand].Path = path;
        fPropagCand[fNbPropagCand].propagState = path->GetFrontSeedState();
        fPropagCand[fNbPropagCand].FrontOrBack = 0;
        fNbPropagCand++;
        if ( ND::tpcDebug().MatchAndMerge(DB_VERBOSE))
          std::cout<<"    PropagCand, Front state, path id: "<<path->GetId()<<std::endl;
      }
      if (SelectedClu == 1 || SelectedClu == 2){
        fPropagCand[fNbPropagCand].Pattern = pattern;
        fPropagCand[fNbPropagCand].Path = path;
        fPropagCand[fNbPropagCand].propagState = path->GetBackSeedState();
        fPropagCand[fNbPropagCand].FrontOrBack = 1;
        fNbPropagCand++;
        if ( ND::tpcDebug().MatchAndMerge(DB_VERBOSE))
          std::cout<<"    PropagCand, Back state, path id: "<<path->GetId()<<std::endl;
      }
    }
  }
}


//*****************************************************************************
bool ND::TTPCSeed2CluMatchMerge::FindMatchCandidates(ND::THandle<ND::TTPCPattern> pattern, int zone){
  if ( ND::tpcDebug().MatchAndMerge(DB_VERBOSE))
    std::cout<<" ---- FindMatchCandidates"<<std::endl;
  if ( ! pattern->IsUsable())
    return false;
  // For each pattern find the paths and junctions at the edge of the MM. Fill simple containers without pattern in it.
  bool FoundCandidate = false;
  for (ND::TReconObjectContainer::iterator constit = pattern->GetConstituents()->begin(); constit != pattern->GetConstituents()->end(); constit++) {
    // Here prevent creation of too many propag candidates.
    if ( !(fNbMatchCand < MAXNBMATCHCANDIDATE))
      break;
    ND::THandle<ND::TTPCPath> path = (*constit);
    if (!path) continue;
    if (!path->GetIsXPath()){

      int SelectedClu;
      GetPathEndsAtZoneEdge(path, SelectedClu, zone);

      if (SelectedClu > -1){
        if ( ND::tpcDebug().MatchAndMerge(DB_VERBOSE))
          std::cout<<"    Selected cluster index: "<<SelectedClu<<std::endl;

        FoundCandidate = true;
        fMatchCand[fNbMatchCand].PathOrJunction = path;
        if ( ND::tpcDebug().MatchAndMerge(DB_VERBOSE))
          std::cout<<"    Selected path id: "<<path->GetId()<<std::endl;
        unsigned int NbClu = 0;
        if (SelectedClu == 0 || SelectedClu == 2){
          fMatchCand[fNbMatchCand].ClusterToMatch[NbClu] = *(path->GetHits()->begin());
          if ( ND::tpcDebug().MatchAndMerge(DB_VERBOSE)){
            TVector3 Pos = fMatchCand[fNbMatchCand].ClusterToMatch[NbClu]->GetPosition();
            std::cout<<"     Matching candidate: first cluster at ("<<Pos.X()<<", "<<Pos.Y()<<", "<<Pos.Z()<<")"<<std::endl;
          }
          NbClu++;
        } else if (SelectedClu == 1 || SelectedClu == 2){
          fMatchCand[fNbMatchCand].ClusterToMatch[NbClu] = *(path->GetHits()->rbegin());
          if ( ND::tpcDebug().MatchAndMerge(DB_VERBOSE)){
            TVector3 Pos = fMatchCand[fNbMatchCand].ClusterToMatch[NbClu]->GetPosition();
            std::cout<<"     Matching candidate: last cluster at ("<<Pos.X()<<", "<<Pos.Y()<<", "<<Pos.Z()<<")"<<std::endl;
          }
          NbClu++;
        }
        fMatchCand[fNbMatchCand].NbClusterToMatch = NbClu;
        // These residuals define the minimum distance between propagated seed and
        // a cluster needed to consider that we have a match
        fMatchCand[fNbMatchCand].Residuals.SetX(fMinDistForMatchP2P);
        fMatchCand[fNbMatchCand].Residuals.SetY(fMinDistForMatchP2P);
        fMatchCand[fNbMatchCand].Residuals.SetZ(fMinDistForMatchP2P);
        fNbMatchCand++;
      }
    }
  }
  for (ND::TReconObjectContainer::iterator constit = pattern->GetConstituents()->begin(); constit != pattern->GetConstituents()->end(); constit++) {
    // Here prevent creation of too many propag candidates.
    if ( !(fNbMatchCand < MAXNBMATCHCANDIDATE))
      break;
    ND::THandle<ND::TTPCJunction> junction = (*constit);
    if (junction) {
      bool SaveJunction = false;
      ND::THitSelection::const_iterator tHit = junction->GetHits()->begin();
      ND::THandle<ND::TTPCHitPad> hPad = *tHit;
      for (; tHit != junction->GetHits()->end(); tHit++){
        hPad = *tHit;
        if (!hPad){
          // TODO: Proper exception
          throw;
        }
        if( HitSelected(hPad,zone)){
          SaveJunction = true;
          break;
        }
      }
      if (SaveJunction){
        FoundCandidate = true;
        // These residuals define the minimum distance between propagated seed and
        // a cluster needed to consider that we have a match
        fMatchCand[fNbMatchCand].Residuals.SetX(fMinDistForMatchP2J);
        fMatchCand[fNbMatchCand].Residuals.SetY(fMinDistForMatchP2J);
        fMatchCand[fNbMatchCand].Residuals.SetZ(fMinDistForMatchP2J);
        fMatchCand[fNbMatchCand].PathOrJunction = junction;
        // Calculate some average position that will be matched to the seeds later on.
        // Do not do a charge average position to avoid a highly ionizing track from
        // moving the position away from the MIPs.
        TVector3 MeanPos(0.0,0.0,0.0);
        for (ND::THitSelection::const_iterator tHit = junction->GetHits()->begin(); tHit != junction->GetHits()->end(); tHit++){
          // We must be careful to take into account the T0
          MeanPos.SetX(MeanPos.X() + junction->GetCalibX((*tHit)));
          MeanPos.SetY(MeanPos.Y() + (*tHit)->GetPosition().Y());
          MeanPos.SetZ(MeanPos.Z() + (*tHit)->GetPosition().Z());
        }
        MeanPos *= 1./(double)junction->GetHits()->size();
        fMatchCand[fNbMatchCand].JuncPosToMatch = MeanPos;
        fNbMatchCand++;

      }
    }
  }
  return FoundCandidate;
}



//*****************************************************************************
void ND::TTPCSeed2CluMatchMerge::GetPathEndsAtZoneEdge(ND::THandle<ND::TTPCPath> path, int &FirstLastAll, int &Zone){
  // FirstLastAll = 0 means first cluster
  // FirstLastAll = 1 means last cluster
  // FirstLastAll = 2 means both first and last cluster
  FirstLastAll = -1;

  ND::THandle<ND::TTPCHVCluster> HVClu[2];
  // Check only the cluster if this side of the track
  // is NOT connected to a junction.
  if ( !path->IsFrontConnected())
    HVClu[0] = *(path->GetHits()->begin());
  if ( !path->IsBackConnected())
    HVClu[1] = *(path->GetHits()->rbegin());

  for (int fl = 0; fl < 2; fl++){
    if ( !HVClu[fl])
      continue;
    if ( !IsClusterUsable(HVClu[fl]))
      continue;
    for (ND::THitSelection::const_iterator tHit = HVClu[fl]->GetHits().begin(); tHit != HVClu[fl]->GetHits().end(); tHit++){
      ND::THandle<ND::TTPCHitPad> hPad = *tHit;
      if (!hPad){
        // TODO: Proper exception
        throw;
      }
      if( HitSelected(hPad,Zone)){
        if (FirstLastAll < 0)
          FirstLastAll = fl;
        else
          FirstLastAll = 2;
        break;
      }
    }
  }
}




//*****************************************************************************
void ND::TTPCSeed2CluMatchMerge::FindMatches(int matchId, ND::THandle<ND::TTPCPattern> MatchCandPat){
  // Loop over the MatchCandidates
  //   Loop over the fPathToPropag
  //     Propagate fPathToPropag to the cluster or pads
  //      if distance less than cutoff (10cm ?) then compare with previous residuals
  //      if residuals smaller than saved one, then save the new one.

  TTPCRecPackUtils::SavedParam RPParam = TTPCRecPackUtils::InitForPropagationInTPC();

  // We have to force the sense of the propagation here for
  // for curving back tracks in particular.

  // If hit to match is HVCluster => easy
  if ( fMatchCand[matchId].NbClusterToMatch ) {
    if ( ND::tpcDebug().MatchAndMerge(DB_VERBOSE))
      std::cout<<" ---- FindMatches to HV cluster from a path"<<std::endl;
    for ( unsigned int nclu = 0; nclu < fMatchCand[matchId].NbClusterToMatch; nclu++){
      if ( ND::tpcDebug().MatchAndMerge(DB_VERBOSE)){
        ND::THandle<ND::TTPCPath> tmpPath = fMatchCand[matchId].PathOrJunction;
        if (!tmpPath)
          continue;
        std::cout<<"    -> MatchCand, path Id "<<tmpPath->GetId()<<", cluster "<<nclu<<std::endl;
      }
      ND::THandle<ND::TTPCHVCluster> Cluster = fMatchCand[matchId].ClusterToMatch[nclu];
      TVector3 clusterPosition = Cluster->GetCalibratedPosition();
      double bestDistance = fMatchCand[matchId].Residuals.Mag();
      for (int ppcd = 0; ppcd < fNbPropagCand; ppcd++){
        if ( ! fPropagCand[ppcd].Pattern->IsUsable())
          continue;
        // Don't match a pattern with itself.
        // Important when the pattern is on top of a MM horizontal gap.
        if ( MatchCandPat == fPropagCand[ppcd].Pattern)
          continue;
        if ( ND::tpcDebug().MatchAndMerge(DB_VERBOSE)){
          std::cout<<"    -> PropagCand, path Id "<<fPropagCand[ppcd].Path->GetId()<<", propagate the ";
          if (fPropagCand[ppcd].FrontOrBack)
            std::cout<<"back state." <<std::endl;
          else
            std::cout<<"front state." <<std::endl;
        }
        int propSense = 1;
        if (!fPropagCand[ppcd].FrontOrBack){
          // Back state => Go forward
          // Front state => Go backward
          propSense = -1;
        }
        ND::rpman("TREx").model_svc().model().intersector().set_length_sign(propSense);
        State propagState = fPropagCand[ppcd].propagState;
        if (!TTPCRecPackUtils::PropagateToHVCluster(Cluster, propagState))
          continue;
        TVector3 propagPosition(propagState.vector()[0], propagState.vector()[1], propagState.vector()[2]);
        // Check the distance
        double Distance = (propagPosition - clusterPosition).Mag();
        // Check that the propagate state sense matches the crude sense given by the matched cluster and the next one in the track.
        // Very helpful to prevent matches of colinear tracks.
        int targetSense = TTPCUtils::SenseFromTwoClusters(fMatchCand[matchId].PathOrJunction, Cluster);
        bool senseOk = true;
        if (  Cluster->IsVertical()  && (propSense*propagState.vector()[5]*targetSense < 0))
          senseOk = false;
        else if ((!Cluster->IsVertical()) && (propSense*propagState.vector()[4]*targetSense < 0))
          senseOk = false;
        
        if (Distance < bestDistance && senseOk){
          if ( ND::tpcDebug().MatchAndMerge(DB_VERBOSE))
            std::cout<<"       MATCH FOUND ==> Distance = "<<Distance<<std::endl;
          fMatchCand[matchId].MatchedPathId = ppcd;
          fMatchCand[matchId].MatchedClusterIdx = nclu;
          fMatchCand[matchId].Residuals = (propagPosition - clusterPosition);
          bestDistance = Distance;
          
        } else {
          if ( ND::tpcDebug().MatchAndMerge(DB_VERBOSE))
            std::cout<<"       MATCH FAILED ==> Distance = "<<Distance<<std::endl;
        }
      }
    }
  } else {
    if ( ND::tpcDebug().MatchAndMerge(DB_VERBOSE)){
      std::cout<<" ---- FindMatches to hit pads from a junction"<<std::endl;
      ND::THandle<ND::TTPCJunction> Junc = fMatchCand[matchId].PathOrJunction;
      if (Junc)
        std::cout<<"    -> MatchCand, junction Id "<<Junc->GetId()<<std::endl;
    }
    double bestDistance = fMatchCand[matchId].Residuals.Mag();
    for (int ppcd = 0; ppcd < fNbPropagCand; ppcd++){
      if ( ! fPropagCand[ppcd].Pattern->IsUsable())
        continue;
      // Don't match a pattern with itself.
      // Important when the pattern is on top of a MM horizontal gap.
      if ( MatchCandPat == fPropagCand[ppcd].Pattern)
        continue;
      if ( ND::tpcDebug().MatchAndMerge(DB_VERBOSE)){
        std::cout<<"    -> PropagCand, path Id "<<fPropagCand[ppcd].Path->GetId()<<", propagate the ";
        if (fPropagCand[ppcd].FrontOrBack)
          std::cout<<"back state." <<std::endl;
        else
          std::cout<<"front state." <<std::endl;
      }
      if (fPropagCand[ppcd].FrontOrBack){  // Back state => Go forward
        ND::rpman("TREx").model_svc().model().intersector().set_length_sign(1);
      } else {                             // Front state => Go backward
        ND::rpman("TREx").model_svc().model().intersector().set_length_sign(-1);
      }
      State propagState = fPropagCand[ppcd].propagState;
      // First propagate to a vertical surface and check the local angle.
      TVector3 projNorm(0.,0.,1.);
      State newState;
      bool RPsucceeded = ND::tman("TREx").PropagateToSurface(propagState,fMatchCand[matchId].JuncPosToMatch,projNorm,newState);
      // Default true so that we redo if RPsucceeded is false.
      bool useHoriSurf = true;
      if (RPsucceeded){
        // The angle is vertical enough. we need a horizontal surface.
        useHoriSurf = (newState.vector()[2] > newState.vector()[1]);
      } 
      // Recalculate the proper position of the propagation using horizontal surface if needed.
      if (useHoriSurf){
        projNorm.SetY(1.);
        projNorm.SetZ(0.);
        RPsucceeded = ND::tman("TREx").PropagateToSurface(propagState,fMatchCand[matchId].JuncPosToMatch,projNorm,newState);
      }
      if (!RPsucceeded)
        continue;
      TVector3 propagPosition(newState.vector()[0], newState.vector()[1], newState.vector()[2]);
      double Distance = (propagPosition - fMatchCand[matchId].JuncPosToMatch).Mag();
      if (Distance < bestDistance){
        if ( ND::tpcDebug().MatchAndMerge(DB_VERBOSE))
          std::cout<<"       MATCH FOUND ==> Distance = "<<Distance<<std::endl;
        fMatchCand[matchId].MatchedPathId = ppcd;
        fMatchCand[matchId].Residuals = (propagPosition - fMatchCand[matchId].JuncPosToMatch);
        bestDistance = Distance;
      } else {
        if ( ND::tpcDebug().MatchAndMerge(DB_VERBOSE))
          std::cout<<"       MATCH FAILED ==> Distance = "<<Distance<<std::endl;
      }
    }
  }
  TTPCRecPackUtils::ResetAfterPropagationInTPC(RPParam);

}


//*****************************************************************************
ND::THandle<ND::TTPCPattern> ND::TTPCSeed2CluMatchMerge::MergePatterns(ND::THandle<ND::TTPCPattern> origPattern, int matchId){
  if ( ND::tpcDebug().MatchAndMerge(DB_VERBOSE))
    std::cout<<" ---- MergePatterns"<<std::endl;
  ND::THandle<ND::TTPCPattern> newPattern;

  int ppcd = fMatchCand[matchId].MatchedPathId;
  ND::THandle<ND::TTPCPath> propPath = fPropagCand[ppcd].Path;
  ND::THandle<ND::TTPCPath> matchPath = fMatchCand[matchId].PathOrJunction;
  ND::THandle<ND::TTPCJunction> matchJunction = fMatchCand[matchId].PathOrJunction;
  ///////////////// path + path
  if ( matchPath) {

    // ====> 1) Paths merging
    ND::THandle<ND::TTPCPath> newPath = TTPCUtils::MergePaths( matchPath, propPath);

    // Redo the seeding from scratch
    fNewSeeding->FindSeed(newPath);

    // TODO: Should we simply recalculate the seed directly here ????
    // Then we wouldn't have to do it later from TRExReco
    
    // Make sure we don't reuse these patterns for another merging.
    origPattern->SetUsable(false);
    fPropagCand[ppcd].Pattern->SetUsable(false);

    // ====> 2) Pattern merging
    if ( ND::tpcDebug().MatchAndMerge(DB_INFO))
      std::cout<<"   Merging MatchCand path "<<matchPath->GetId()<<" and PropagCand path "<<propPath->GetId()<<" into the new path "<<newPath->GetId()<<std::endl;

    //   IF two single tracks, DONE HERE !
    if ( matchPath->GetNJunctionIds() == 0 && propPath->GetNJunctionIds() == 0){
      if (fFakeJunctionForStudies && fAlgo == kMMVertGapMerge){
        newPattern = ND::THandle<ND::TTPCPattern>( new ND::TTPCPattern() );
        ND::THandle<ND::TTPCHVCluster> RefClu;
        if (fPropagCand[ppcd].FrontOrBack)
          RefClu = *(propPath->GetHits()->rbegin());
        else
          RefClu = *(propPath->GetHits()->begin());
        ND::THandle<ND::THit> FirstHit = *(RefClu->GetHits().begin());
        TLorentzVector NewPos(RefClu->GetCalibPosition(),RefClu->GetTime());
        ND::THandle<ND::TTPCJunction> newJunction = ND::THandle<ND::TTPCJunction>( new ND::TTPCJunction(NewPos) );
        ND::THitSelection *hitJunction = new ND::THitSelection();
        hitJunction->AddHit(FirstHit);
        newJunction->AddHits(hitJunction);
        newJunction->AddConstituent(matchPath);
        newJunction->AddConstituent(propPath);
        newPattern->AddJunction(newJunction);
        newPattern->SetId(ND::tpcCalibration().GetPatternId());
        newPattern->InitialSetup();

        // Make sure we don't reuse these patterns for another merging.
        origPattern->SetUsable(false);
        fPropagCand[ppcd].Pattern->SetUsable(false);
      } else {
        if ( ND::tpcDebug().MatchAndMerge(DB_VERBOSE))
          std::cout<<"   -> Simple track-track merging"<<std::endl;
        newPattern = ND::THandle<ND::TTPCPattern>( new ND::TTPCPattern(newPath) );
        newPattern->SetId(ND::tpcCalibration().GetPatternId());
        newPattern->InitialSetup();
      }

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

      if ( ND::tpcDebug().MatchAndMerge(DB_VERBOSE)){
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
        MigrateJunctions(oldPattern[nC], oldPath[nC]->GetId(), newPath, newPattern);
      }
      newPattern->SetId(ND::tpcCalibration().GetPatternId());
      newPattern->InitialSetup();
    }

    // if the match candidate had more than one cluster to match, then get rid of the one we used.
    if (fMatchCand[matchId].NbClusterToMatch > 1){
      fMatchCand[matchId].NbClusterToMatch = 1;
      if (fMatchCand[matchId].MatchedClusterIdx == 0)
        fMatchCand[matchId].ClusterToMatch[0] = fMatchCand[matchId].ClusterToMatch[1];
    } else {
      fMatchCand[matchId].NbClusterToMatch = 0;
      fMatchCand[matchId].ClusterToMatch[0] = ND::THandle<ND::TTPCHVCluster>();
    }
    // properly reset the match variables
    fMatchCand[matchId].ClusterToMatch[1] = ND::THandle<ND::TTPCHVCluster>();
    fMatchCand[matchId].Residuals = TVector3(fMinDistForMatchP2P, fMinDistForMatchP2P, fMinDistForMatchP2P);
    fMatchCand[matchId].MatchedPathId = -1;
    fMatchCand[matchId].MatchedClusterIdx = -1;

    return newPattern;
  }
  ///////////////// path + junction
  else if (matchJunction){
    fPropagCand[ppcd].Pattern->SetUsable(false);
    if ( ND::tpcDebug().MatchAndMerge(DB_INFO))
      std::cout<<"   Merging MatchCand junction "<<matchJunction->GetId()<<" and PropagCand path "<<propPath->GetId()<<std::endl;
    // Just add the propagated path to the matched junction.
    for (ND::TReconObjectContainer::iterator constit = origPattern->GetConstituents()->begin(); constit != origPattern->GetConstituents()->end(); constit++) {
      ND::THandle<ND::TTPCJunction> origJunction = (*constit);
      if (!origJunction)
        continue;
      if ( origJunction == matchJunction){
        origJunction->AddConstituent(propPath);
        // Do it manually because we won't call AddJunction with origJunction.
        origPattern->AddConstituent(propPath);
        break;
      }
    }
    // If the propagated path comes from a single clean track pattern, then we are done.
    if ( propPath->GetNJunctionIds() == 0){
      origPattern->InitialSetup();
      return origPattern;
    }
    // Migrate the junctions from the propagated pattern to the matched pattern to extend the latter.
    else {
      for (ND::TReconObjectContainer::iterator constit = fPropagCand[ppcd].Pattern->GetConstituents()->begin(); constit != fPropagCand[ppcd].Pattern->GetConstituents()->end(); constit++) {
        ND::THandle<ND::TTPCJunction> propJunction = (*constit);
        if (!propJunction)
          continue;
        origPattern->AddJunction(propJunction);
      }
    }
    origPattern->InitialSetup();
    return origPattern;
  }

  std::cerr<<" This is a very bad day for you. The matched object is neither a TTPCPath, nor a TTPCJunction. Let's crash !"<<std::endl;
  throw;
}


//*****************************************************************************
void ND::TTPCSeed2CluMatchMerge::MigrateJunctions(ND::THandle<ND::TTPCPattern> PatternA, unsigned int PathIdA, ND::THandle<ND::TTPCPath> PathAB, ND::THandle<ND::TTPCPattern> PatternB){
  // Recreate the junctions attached to the merged paths, making sure
  // we replace the matched paths by the new path in the transfer of the paths
  // from the old to the new junctions.
  for (ND::TReconObjectContainer::iterator constit = PatternA->GetConstituents()->begin(); constit != PatternA->GetConstituents()->end(); constit++) {
    // Find the junction that is connected to the merged path in the old pattern.
    ND::THandle<ND::TTPCJunction> JunctionA = (*constit);
    if (!JunctionA)
      continue;
    bool HasMergedPath = false;
    for (ND::TReconObjectContainer::iterator pathtit = JunctionA->GetConstituents()->begin(); pathtit != JunctionA->GetConstituents()->end(); pathtit++) {
      ND::THandle<ND::TTPCPath> tmpPath = (*pathtit);
      if ( tmpPath->GetId() == PathIdA){
        HasMergedPath = true;
        break;
      }
    }
    if (HasMergedPath) {
      // Need JunctionA position
      ND::THandle<ND::TTPCJunction> newJunction( new ND::TTPCJunction(JunctionA->GetPosition()) );
      ND::THandle<ND::THitSelection> oldHits = JunctionA->GetHits();
      ND::THitSelection *newHits = new ND::THitSelection();
      for (ND::THitSelection::const_iterator tmpHit = oldHits->begin(); tmpHit != oldHits->end(); tmpHit++) {
        newHits->push_back(*tmpHit);
      }
      newJunction->AddHits(newHits);
      for (ND::TReconObjectContainer::iterator pathtit = JunctionA->GetConstituents()->begin(); pathtit != JunctionA->GetConstituents()->end(); pathtit++) {
        ND::THandle<ND::TTPCPath> tmpPath = (*pathtit);
        if ( tmpPath->GetId() == PathIdA){
          newJunction->AddConstituent(PathAB);
        } else {
          newJunction->AddConstituent(tmpPath);
        }
      }
      newJunction->SetId(ND::tpcCalibration().GetJunctionId());
      if ( ND::tpcDebug().MatchAndMerge(DB_VERBOSE))
        std::cout<<"        New junction created with Id "<<newJunction->GetId()<<std::endl;
      PatternB->AddJunction(newJunction);
    } else {
      PatternB->AddJunction(JunctionA);
    }
  }
}



