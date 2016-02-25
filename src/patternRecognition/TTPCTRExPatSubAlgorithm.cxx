// eddy
#include "TTPCTRExPatSubAlgorithm.hxx"
#include "TTPCUtils.hxx"

ND::TTPCTRExPatSubAlgorithm::TTPCTRExPatSubAlgorithm(ND::TTPCLayout* layout) {
  // initial values for containing hits, tpc number and event and sub-event numbers
  fHasHits = false;
  fHasValidPaths = false;
  fPrimary = false;
  fTPC = 0;

  // initialise empty vectors
  fPOIs = std::vector<TVector3>();
  fGOIs = std::vector< ND::THandle<ND::TTPCVolGroup> >();
  fTracks = std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >();

  // set up objects for layout, hit group manager, feature finder and path finder
  fLayout = layout;

  fVolGroupMan = new ND::TTPCVolGroupMan(fLayout);
  fVolGroupManHigh = new ND::TTPCVolGroupMan(fLayout);
  fAStar = new ND::TTPCAStar(fLayout);
  fAStarHigh = new ND::TTPCAStar(fLayout);

  fPattern = ND::THandle< ND::TTPCPattern >(); //0;
}
ND::TTPCTRExPatSubAlgorithm::~TTPCTRExPatSubAlgorithm(){
  // delete groups of hits and path finders
  delete fVolGroupMan;
  delete fVolGroupManHigh;
  delete fAStar;
  delete fAStarHigh;
}

void ND::TTPCTRExPatSubAlgorithm::SetUpHits(std::map<long, ND::TTPCUnitVolume*> map, ND::TTPCAStar* aStarCopy){
  // add new pattern recognition cells if the charge cut is met
  fHitMap = std::map<long, ND::TTPCUnitVolume*>(map);
  fHasHits = fHitMap.size() > 0;

  // add map of hits to group manager object
  fVolGroupMan->AddPrimaryHits(fHitMap);
  // add map of hits to path finder
  aStarCopy = 0;
  if(aStarCopy){
    fAStar->AddHits(aStarCopy, fHitMap);
  }
  else{
    fAStar->AddHits(fVolGroupMan, fHitMap);
  };
}
void ND::TTPCTRExPatSubAlgorithm::SetUpHitsHigh(std::map<long, ND::TTPCUnitVolume*> map, ND::TTPCAStar* aStarCopy){
  // add new pattern recognition cells if the charge cut is met
  fHitMapHigh = std::map<long, ND::TTPCUnitVolume*>(map);

  // add map of high charge hits to group manager object
  fVolGroupManHigh->AddPrimaryHits(fHitMapHigh);
  // add map of high charge hits to path finder
  aStarCopy = 0;
  if(aStarCopy){
    fAStarHigh->AddHits(aStarCopy, fHitMapHigh);
  }
  else{
    fAStarHigh->AddHits(fVolGroupManHigh, fHitMapHigh);
  };
}
void ND::TTPCTRExPatSubAlgorithm::AbsorbCell(long id, ND::TTPCUnitVolume* cell){
  // if cell doesn't already exist in map, create it; if it does, increment its charge by cell's charge
  std::map<long, ND::TTPCUnitVolume*>::iterator el = fHitMap.find(id);
  if(el == fHitMap.end()){
    fHitMap[id] = cell;
  }
  else{
    fHitMap[id]->AddCharge(cell->GetQ());
    fHitMap[id]->AddHits(cell->GetHits());
    delete cell;
  };
}
void ND::TTPCTRExPatSubAlgorithm::AppendHits(ND::THandle<ND::TTPCVolGroup> hits){
  // loop over all cells in hit group and absorb each one into this
  for(std::map<long, ND::TTPCUnitVolume*>::iterator it = hits->begin(); it != hits->end(); ++it){
    AbsorbCell(it->first, it->second);
  };
}
void ND::TTPCTRExPatSubAlgorithm::ProduceContainers(){
  // ignore if hits don't already exist
  if(!fHasHits) return;

  // start out having no valid paths by default
  fHasValidPaths = false;

  // show regions of pathological behaviour for debugging
  /*ND::THandle<ND::TTPCVolGroup> pathologicalGroup (new ND::TTPCVolGroup(fLayout));
  for(std::map<long, ND::TTPCUnitVolume*>::iterator el = fHitMap.begin(); el != fHitMap.end(); ++el){
    if(el->second->GetPathology()){
    //if(el->second->GetSatASICTagged()){
    //if(el->second->GetFullASICTagged()){
    //if(el->second->GetLowChargeTagged()){
    //if(el->second->GetLateNegativeTagged() || el->second->GetEarlyNegativeTagged()){
      pathologicalGroup->AddHit(el->second);
    };
  };
  fGOIs.push_back(pathologicalGroup);*/

  // find effective edge detection vol group and path finder based on charge cut
  // use hybridVolGroupMan/hybridAStar for edge detection
  ND::TTPCVolGroupMan* hybridVolGroupMan;
  ND::TTPCAStar* hybridAStar;
  if(fLayout->GetUsePatRecPathologyCut() > 0){
    hybridVolGroupMan = fVolGroupManHigh;
    hybridAStar = fAStarHigh;
  }
  else{
    hybridVolGroupMan = fVolGroupMan;
    hybridAStar = fAStar;
  };
  // find effective path finding vol group and path finder based on charge cut
  // use effVolGroupMan/effAStar for path finding and fVolGroupMan/fAStar for hit association
  ND::TTPCVolGroupMan* effVolGroupMan;
  ND::TTPCAStar* effAStar;
  if(fLayout->GetUsePatRecPathologyCut() == 1){
    effVolGroupMan = fVolGroupManHigh;
    effAStar = fAStarHigh;
  }
  else{
    effVolGroupMan = fVolGroupMan;
    effAStar = fAStar;
  };

  // find list of hit groups on edge of paths
  std::vector< ND::THandle<ND::TTPCVolGroup> > edgeGroupsUnfiltered;
  if(fLayout->GetUseAltEdgeDetection()) edgeGroupsUnfiltered = hybridVolGroupMan->GetAllEdgeGroups();
  else edgeGroupsUnfiltered = hybridVolGroupMan->GetEdgeGroups();

  // clear redundant edge groups (i.e. ones which are actually just half way points along existing paths)
  std::vector< ND::THandle<ND::TTPCVolGroup> > edgeGroupsSemiFiltered;
  if(fLayout->GetUsePatRecPathologyCut() > 1) edgeGroupsSemiFiltered = effAStar->ExperimentalClearRedundancies2(effVolGroupMan, edgeGroupsUnfiltered, true);
  else edgeGroupsSemiFiltered = effAStar->ClearRedundancies(effVolGroupMan, edgeGroupsUnfiltered);

  // connect pairs of edge groups to find missing hits
  std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > edgePathsSemiFiltered = effAStar->ConnectGroupsOrdered(effVolGroupMan, edgeGroupsSemiFiltered, false, false);

  // find unmatched hits
  std::map<long, ND::TTPCUnitVolume*> unmatchedVolumes;
  if(fLayout->GetUsePatRecPathologyCut() > 0){
    unmatchedVolumes= std::map<long, ND::TTPCUnitVolume*>(fHitMapHigh);
  }
  else{
    unmatchedVolumes = std::map<long, ND::TTPCUnitVolume*>(fHitMap);
  };
  if(fLayout->GetUseAltHitAssociation() > 0) effAStar->AssociateBestHits(effVolGroupMan, edgePathsSemiFiltered, fLayout->GetAltExtraHitConnectDist());
  for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator edgePathIt = edgePathsSemiFiltered.begin(); edgePathIt != edgePathsSemiFiltered.end(); ++edgePathIt){
    ND::THandle<ND::TTPCOrderedVolGroup> edgePath = *edgePathIt;

    if(fLayout->GetUseAltHitAssociation() == 0) hybridVolGroupMan->BuildGroupFriends(edgePath, ND::TTPCConnection::extraHits);
    ND::THandle<ND::TTPCVolGroup> pathHits = edgePath->GetExtendedHits();

    for(std::map<long, ND::TTPCUnitVolume*>::iterator foundVolIt = pathHits->begin(); foundVolIt != pathHits->end(); ++foundVolIt){
      long id = foundVolIt->first;
      std::map<long, ND::TTPCUnitVolume*>::iterator foundIdLoc = unmatchedVolumes.find(id);
      if(foundIdLoc != unmatchedVolumes.end()) unmatchedVolumes.erase(foundIdLoc);
    };
  };

  // split and group unmatched hits
  ND::TTPCVolGroupMan* secondPassVolume = new ND::TTPCVolGroupMan(fLayout);
  secondPassVolume->AddPrimaryHits(unmatchedVolumes);

  std::vector< ND::THandle<ND::TTPCVolGroup> > extraSubVolumes = secondPassVolume->GetConnectedHits(ND::TTPCConnection::path, ND::TTPCHitGroupings::all, true);
  for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator extraSubVolumeIt = extraSubVolumes.begin(); extraSubVolumeIt != extraSubVolumes.end(); ++extraSubVolumeIt){
    ND::THandle<ND::TTPCVolGroup> extraSubVolume = *extraSubVolumeIt;

    ND::TTPCVolGroupMan* subVolumeMan = new ND::TTPCVolGroupMan(fLayout);
    std::map<long, ND::TTPCUnitVolume*> subVolumeHits = std::map<long, ND::TTPCUnitVolume*>();
    for(std::map<long, ND::TTPCUnitVolume*>::iterator volIt = extraSubVolume->begin(); volIt != extraSubVolume->end(); ++volIt) subVolumeHits[volIt->first] = volIt->second;

    subVolumeMan->AddPrimaryHits(subVolumeHits);
    std::vector< ND::THandle<ND::TTPCVolGroup> > extraEdgeGroupsSemiFiltered = subVolumeMan->GetEdgeGroups();
    // add unmatched edges to main group
    for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator extraEdgeGroupIt = extraEdgeGroupsSemiFiltered.begin(); extraEdgeGroupIt != extraEdgeGroupsSemiFiltered.end(); ++extraEdgeGroupIt){
      edgeGroupsSemiFiltered.push_back(*extraEdgeGroupIt);
    };

    delete subVolumeMan;
  };
  delete secondPassVolume;

  // clear redundant edge groups (i.e. ones which are actually just half way points along existing paths)
  std::vector< ND::THandle<ND::TTPCVolGroup> > edgeGroups;
  if(fLayout->GetUsePatRecPathologyCut() > 1) edgeGroups = effAStar->ExperimentalClearRedundancies2(effVolGroupMan, edgeGroupsSemiFiltered, true);
  else edgeGroups = effAStar->ClearRedundancies(effVolGroupMan, edgeGroupsSemiFiltered);

  // add to points of interest for debugging
  //for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator edgeGroup = edgeGroupsUnfiltered.begin(); edgeGroup != edgeGroupsUnfiltered.end(); ++edgeGroup){
  //for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator edgeGroup = edgeGroupsSemiFiltered.begin(); edgeGroup != edgeGroupsSemiFiltered.end(); ++edgeGroup){
  //for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator edgeGroup = edgeGroups.begin(); edgeGroup != edgeGroups.end(); ++edgeGroup){
  //  fPOIs.push_back((*edgeGroup)->GetAveragePosXYZ());
  //  fGOIs.push_back(*edgeGroup);
  //};

  // connect pairs of edge groups
  std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > edgePaths = effAStar->ConnectGroupsOrdered(effVolGroupMan, edgeGroups);

  // define container for true paths before corner detection
  std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > truePaths = std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >();

  // define vertices for later use
  std::vector< ND::THandle<ND::TTPCVolGroup> > vertices = std::vector< ND::THandle<ND::TTPCVolGroup> >();

  // if there are at least three edge groups, look for vertices to connect to them
  if(edgeGroups.size() > 2){
    // connect pairs of edge groups
    std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > edgePaths = effAStar->ConnectGroupsOrdered(effVolGroupMan, edgeGroups, false, true);

    // look for vertices from edge group pairs
    if(fLayout->GetUsePatRecPathologyCut() > 1){
      std::vector< ND::THandle<ND::TTPCVolGroup> > verticesUnfiltered = effVolGroupMan->GetFoci(edgePaths, .1);
      // clean up vertices near to each other
      vertices = effAStar->ExperimentalClearVertexRedundancies(verticesUnfiltered);
    }
    else{
      vertices = effVolGroupMan->GetFoci(edgePaths, .1);
    };

    // make sure number of vertices is two less than number of track ends
    effVolGroupMan->CleanUpVertices(edgeGroups, vertices);
    // make sure vertices contain all the right hits
    fVolGroupMan->BulkGroups(vertices);

    //for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator vertex = vertices.begin(); vertex != vertices.end(); ++vertex){
    //  fPOIs.push_back((*vertex)->GetAveragePosXYZ());
    //  fGOIs.push_back(*vertex);
    //};

    if(vertices.size() > 0){
      // connect vertices to path ends to get paths
      std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > truePathsUnfiltered = effAStar->ConnectVertexGroupsOrdered(effVolGroupMan, vertices, edgeGroups);
      // clear redundant paths
      truePaths = effAStar->ClearVertexConnectionRedundancies(effVolGroupMan, truePathsUnfiltered, vertices);
    }
    // otherwise, return empty handed
    else{
      return;
    };
  }
  // if there are two, connect the two edge groups
  else if (edgeGroups.size() > 1){
    truePaths = edgePaths;
  }
  // otherwise, return empty handed
  else{
    return;
  };

  // ensure that paths don't share hits with each other
  if(fLayout->GetUseAltHitAssociation() >= 1) hybridAStar->AssociateBestHits(fVolGroupMan, truePaths, fLayout->GetAltPathHitConnectDist());
  else fVolGroupMan->BuildAllFriends(truePaths);

  // extend and cluster paths
  for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator truePathIt = truePaths.begin(); truePathIt != truePaths.end(); ++truePathIt){
    fVolGroupMan->ClusterGroupFriends(*truePathIt, false, true);
  };
  // clear empties
  std::vector< ND::THandle<ND::TTPCOrderedVolGroup > >::iterator pathDel = truePaths.begin();
  while(pathDel != truePaths.end()){
    if(!(*pathDel)->size()) truePaths.erase(pathDel);
    else pathDel ++;
  };

  // look for kinks in paths
  std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > brokenPaths = hybridVolGroupMan->BreakPathsAboutKinks(truePaths, true);
  if(fLayout->GetUseAltHitAssociation() >= 1) hybridAStar->AssociateBestHits(effVolGroupMan, brokenPaths, fLayout->GetAltPathHitConnectDist());

  std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > brokenClusteredPaths;
  for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator brokenPathIt = brokenPaths.begin(); brokenPathIt != brokenPaths.end(); ++brokenPathIt){
    ND::THandle<ND::TTPCOrderedVolGroup> brokenPath = *brokenPathIt;

    effVolGroupMan->ClusterGroupFriends(brokenPath, true, true);
    brokenClusteredPaths.push_back(brokenPath);
  };

  // create list of final paths
  std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > finalPaths;
  for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator pathIt = brokenClusteredPaths.begin(); pathIt != brokenClusteredPaths.end(); ++pathIt){
    finalPaths.push_back(*pathIt);
  };

  // temporarily merge x-paths into junctions to make the next two steps easier
  effVolGroupMan->SeparateXPathHits(finalPaths);

  // clean up any and all instances of hits beiing placed in the wrong object
  effVolGroupMan->SeparateEmptyClusters(finalPaths);
  effVolGroupMan->SeparateClusterHits(finalPaths);
  effVolGroupMan->SeparateAnomHits(finalPaths);
  effVolGroupMan->SeparateJunctionHits(finalPaths);

  // enforce path ordering
  effVolGroupMan->EnforceOrdering(finalPaths);

  // add all paths that make sense and store if none are found
  fTracks = effVolGroupMan->SanityFilter(finalPaths);

  // associate any remaining unused hits to junctions where possible
  ND::THandle<ND::TTPCVolGroup> unusedHits = effVolGroupMan->GetUnusedHits(fTracks);
  effVolGroupMan->AssociateUnusedWithJunctions(unusedHits, fTracks);

  // break junctions with a large extent in x
  if(!fLayout->GetUseAltEdgeDetection()) effVolGroupMan->BreakLongJunctions(fTracks);

  // reset vertex status based on how many paths a junction is connected to
  effVolGroupMan->ResetVertexStatuses(fTracks);

  // if using a charge cut, associate low charge hits with paths and junctions as a final step
  if(fLayout->GetUsePatRecPathologyCut() > 0){
    if(fLayout->GetUseAltHitAssociation() >= 1){
      fAStar->AssociateBestHits(fVolGroupMan, fTracks, fLayout->GetAltPathHitConnectDist());
    }
    else{
      fVolGroupMan->BuildAllFriends(fTracks);
    };

    for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator trackIt = fTracks.begin(); trackIt != fTracks.end(); ++trackIt){
      ND::THandle<ND::TTPCOrderedVolGroup> track = *trackIt;
      fVolGroupMan->ClusterGroupFriends(track, true, true, true);
    };

    effVolGroupMan->AssociateUnusedWithJunctions(unusedHits, fTracks);
  };

  // get junctions from paths for diagnostics
  //std::vector< ND::THandle<ND::TTPCVolGroup> > pathVertices = fVolGroupMan->GetJunctionsFromPaths(fTracks);
  //for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator pathVertexIt = pathVertices.begin(); pathVertexIt != pathVertices.end(); ++pathVertexIt){
  //  ND::THandle<ND::TTPCVolGroup> pathVertex = *pathVertexIt;
  //  fPOIs.push_back(pathVertex->GetAveragePosXYZ());
  //  fGOIs.push_back(pathVertex);
  //};

  fHasValidPaths = fTracks.size() > 0;
}

ND::THandle<ND::TTPCPattern> ND::TTPCTRExPatSubAlgorithm::GetPattern(){
  if(!fPattern){
    ND280Warn("Pattern not set for TTPCTRExPatSubAlgorithm");
    //TF: an empty handle on the pattern, don't know whether this is the best idea
    return ND::THandle< ND::TTPCPattern >(); //0;
  };
  return fPattern;
}
void ND::TTPCTRExPatSubAlgorithm::ProducePattern(ND::THitSelection* used){
  if(fPattern){
    ND280Warn("Trying to create same pattern twice for TTPCTRExPatSubAlgorithm");
    return;
  };

  // final paths and junctions
  std::vector< ND::THandle<ND::TTPCPath> > paths = std::vector< ND::THandle<ND::TTPCPath> >();
  std::vector< ND::THandle<ND::TTPCJunction> > junctions = std::vector< ND::THandle<ND::TTPCJunction> >();

  // temporary vectors for holding candidate paths and junctions, and parallel vector of groups associated with junctions
  std::vector< ND::THandle<ND::TTPCPath> > tempPaths = std::vector< ND::THandle<ND::TTPCPath> >();
  std::vector< ND::THandle<ND::TTPCJunction> > tempJunctions = std::vector< ND::THandle<ND::TTPCJunction> >();
  std::vector< ND::THandle<ND::TTPCVolGroup> > tempJunctionGroups = std::vector< ND::THandle<ND::TTPCVolGroup> >();

  for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator trackIt = fTracks.begin(); trackIt != fTracks.end(); ++trackIt){
    ND::THandle<ND::TTPCOrderedVolGroup> track = *trackIt;

    // ignore if empty
    ND::THandle<ND::THitSelection> clusters = track->GetClusters();
    if(!clusters->size()) continue;

    // create new path
    ND::THandle<ND::TTPCPath> path( new ND::TTPCPath() );
    tempPaths.push_back(path);

    // otherwise add all hits to path
    ND::THitSelection* copyClusters = new ND::THitSelection(*clusters);
    path->AddHits(copyClusters);

    // also save if delta criteria met
    if(track->GetDeltaCriteriaMet()){
      path->SetPID(ND::TReconPID::kEM, 1.);
    };
    // and save if the path is in the x direction
    if(track->GetIsXPath()){
      path->SetIsXPath(true);
    };

    // save groups for junctions
    std::vector< ND::THandle<ND::TTPCVolGroup> > junctionGroups = std::vector< ND::THandle<ND::TTPCVolGroup> >();
    if(track->HasBackHits()) junctionGroups.push_back(track->GetBackHits());
    if(track->HasFrontHits()) junctionGroups.push_back(track->GetFrontHits());
    // determine if junction candidate has already been found and create a group for it if it hasn't
    for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator junctionGroupIt = junctionGroups.begin(); junctionGroupIt != junctionGroups.end(); ++junctionGroupIt){
      ND::THandle<ND::TTPCVolGroup> junctionGroup = *junctionGroupIt;
      bool found = false;

      int iMax = tempJunctions.size();
      for(int i=0; i<iMax; i++){
        ND::THandle<ND::TTPCJunction> tempJunction = tempJunctions[i];
        ND::THandle<ND::TTPCVolGroup> tempJunctionGroup = tempJunctionGroups[i];

        // if the junction has the same id as the temporary, add this path to it
        if(junctionGroup->GetID() == tempJunctionGroup->GetID()){
          // if already found, check this path isn't already added and add this path
          bool pathInJunction = false;
          for(ND::TReconObjectContainer::iterator junctionConstIt = tempJunction->GetConstituents()->begin(); junctionConstIt != tempJunction->GetConstituents()->end(); ++junctionConstIt){
            if(path == *junctionConstIt){
              pathInJunction = true;
              break;
            };
          };
          if(!pathInJunction) tempJunction->AddConstituent(path);

          found = true;
        };
      };
      // if not alread found, create a new one
      if(!found){
        ND::THitSelection* junctionHitSelection = junctionGroup->GetHitSelection();
        TVector3 tempVect3 = junctionGroup->GetAverageVol()->GetPos();
        double tempTime = junctionGroup->GetAverageVol()->GetTime();
        // Get a guess X position based on a default T0 of 0.0 just to have a useful X position downstream.
        double tempX = TTPCUtils::GetXWithT0(tempTime, 0.0, (*(junctionHitSelection->begin())));
        tempVect3.SetX(tempX);
        TLorentzVector tempPosition(tempVect3, tempTime);

        ND::THandle<ND::TTPCJunction>  junct( new ND::TTPCJunction(tempPosition) );
        junct->AddHits(junctionHitSelection);
        junct->AddConstituent(path);

        // make sure these two vectors remain parallel
        tempJunctions.push_back(junct);
        tempJunctionGroups.push_back(junctionGroup);
      };
    };
  };

  // add paths that meet criteria
  for(std::vector< ND::THandle<ND::TTPCPath> >::iterator path = tempPaths.begin(); path != tempPaths.end(); ++path){
    (*path)->SetId(ND::tpcCalibration().GetPathId());
    paths.push_back(*path);
  };

  // add junctions that meet criteria
  for(std::vector< ND::THandle<ND::TTPCJunction> >::iterator junction = tempJunctions.begin(); junction != tempJunctions.end(); ++junction){
    if((*junction)->GetNPaths() > 1){
      (*junction)->SetId(ND::tpcCalibration().GetJunctionId());
      junctions.push_back(*junction);
    }
  };

  // add junctions only.  paths will be extracted from the junctions
  if(junctions.size() == 0){
    if(paths.size() == 1){
      ND::THandle<ND::TTPCPath> thepath = *(paths.begin());
      fPattern = ND::THandle<ND::TTPCPattern>( new ND::TTPCPattern(thepath) );
      // One clean single track. Special treatment.
    }
    else{
      // something went wrong - probably an unreconstructable event.  save a null pointer so the main algorithm knows not to add it to list of patterns.
      if(ND::tpcDebug().PatternRecognition(DB_ERROR)) std::cout << "WARNING: no pattern reconstructed " << std::endl;
      fPattern = ND::THandle<ND::TTPCPattern>();
    };
  } else {
    fPattern = ND::THandle<ND::TTPCPattern>( new ND::TTPCPattern() );
    for(std::vector< ND::THandle<ND::TTPCJunction> >::iterator junction = junctions.begin(); junction != junctions.end(); ++junction)
      fPattern->AddJunction(*junction);
  };
  // fill used and unused hits
  FillUsedHits(used);

  // Call the setup method to define things like the detector bit.
  if(fPattern){
    fPattern->SetId(ND::tpcCalibration().GetPatternId());
    fPattern->InitialSetup();
  }
}

ND::THandle<ND::THitSelection> ND::TTPCTRExPatSubAlgorithm::GetHits(){
  return fVolGroupMan->GetHits();
}
ND::THandle<ND::THitSelection> ND::TTPCTRExPatSubAlgorithm::GetHits(ND::THandle<ND::TTPCOrderedVolGroup> path){
  return path->GetClusters();
}

std::vector< ND::THandle<ND::TTPCVolGroup> > ND::TTPCTRExPatSubAlgorithm::GetRegions(){
  if(!fHasHits) return std::vector< ND::THandle<ND::TTPCVolGroup> >();
  // get list of hit groups for sub-events
  std::vector< ND::THandle<ND::TTPCVolGroup> > connectedHits = fVolGroupMan->GetConnectedHits(ND::TTPCConnection::path);

  // for debugging, set as this group's group of interest list
  //fGOIs = connectedHits;

  // return it
  return connectedHits;
}

void ND::TTPCTRExPatSubAlgorithm::FillUsedHits(ND::THitSelection* used){
  // if pattern hasn't been correctly built, hits can't possibly be used
  if(!fPattern) return;

  // matched (initially empty) volumes
  std::map<long, ND::TTPCUnitVolume*> matchedVolumes = std::map<long, ND::TTPCUnitVolume*>();

  // find those contained in any paths or junctions
  for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator trackIt = fTracks.begin(); trackIt != fTracks.end(); ++trackIt){
    ND::THandle<ND::TTPCOrderedVolGroup> track = *trackIt;
    // path hits
    for(std::vector<ND::TTPCPathVolume*>::iterator pathIt = track->begin(); pathIt != track->end(); ++pathIt){
      ND::TTPCPathVolume* path = *pathIt;

      std::vector<ND::TTPCUnitVolume*> pathClusterVols = path->GetExtendedCell(0);
      for(std::vector<ND::TTPCUnitVolume*>::iterator volIt = pathClusterVols.begin(); volIt != pathClusterVols.end(); ++volIt){
        std::map<long, ND::TTPCUnitVolume*>::iterator foundVol = fHitMap.find((*volIt)->GetID());
        matchedVolumes[foundVol->first] = foundVol->second;
      };
    };
  };

  // junction hits
  std::vector< ND::THandle<ND::TTPCVolGroup> > pathVertices = fVolGroupMan->GetJunctionsFromPaths(fTracks);
  for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator pathVertexIt = pathVertices.begin(); pathVertexIt != pathVertices.end(); ++pathVertexIt){
    ND::THandle<ND::TTPCVolGroup> pathVertex = *pathVertexIt;

    for(std::map<long, ND::TTPCUnitVolume*>::iterator volIt = pathVertex->begin(); volIt != pathVertex->end(); ++volIt){
      std::map<long, ND::TTPCUnitVolume*>::iterator foundVol = fHitMap.find(volIt->first);
      matchedVolumes[foundVol->first] = foundVol->second;
    };
  };

  // add used hits
  for(std::map<long, ND::TTPCUnitVolume*>::iterator volEl = matchedVolumes.begin(); volEl != matchedVolumes.end(); ++volEl){
    ND::TTPCUnitVolume* vol = volEl->second;
    if(!vol) continue;

    std::vector< ND::THandle<ND::TTPCHitPad> > hits = vol->GetHits();
    for(std::vector< ND::THandle<ND::TTPCHitPad> >::iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt){
      ND::THandle<ND::TTPCHitPad> hitPad = *hitIt;
      if(hitPad) used->AddHit(hitPad);
    };
  };
}

std::vector< ND::THandle<ND::TTPCVolGroup> > ND::TTPCTRExPatSubAlgorithm::GetTrackExtendedHits(){
  std::vector< ND::THandle<ND::TTPCVolGroup> > outGroups;

  for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator trackIt = fTracks.begin(); trackIt != fTracks.end(); ++trackIt){
    ND::THandle<ND::TTPCOrderedVolGroup> track = *trackIt;
    outGroups.push_back(track->GetExtendedHits());
  };

  return outGroups;
}
