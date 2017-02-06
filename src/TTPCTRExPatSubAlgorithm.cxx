// eddy
#include "TTPCTRExPatSubAlgorithm.hxx"

trex::TTPCTRExPatSubAlgorithm::TTPCTRExPatSubAlgorithm(trex::TTPCLayout* layout) {
  // initial values for containing hits, tpc number and event and sub-event numbers
  fHasHits = false;
  fHasValidPaths = false;
  fPrimary = false;
  fTPC = 0;
  std::cout<<"1"<<std::endl;
  // set up objects for layout, hit group manager, feature finder and path finder
  fLayout = layout;
  std::cout<<"2"<<std::endl;
  fVolGroupMan = new trex::TTPCVolGroupMan(fLayout);
  std::cout<<"2a"<<std::endl;
  fAStar = new trex::TTPCAStar(fLayout);  
  std::cout<<"3"<<std::endl;
  //  fPattern = new trex::TTPCPattern;
}
trex::TTPCTRExPatSubAlgorithm::~TTPCTRExPatSubAlgorithm(){
  // delete groups of hits and path finders
  if(fVolGroupMan) delete fVolGroupMan;
  if(fAStar) delete fAStar;
  //MDH
  //Nobody had better use the pattern after we destroy the subalgorithm...
  //This is an output object that needs to be thoroughly re-engineered anyway
  //delete fPattern;
}

trex::TTPCTRExPatSubAlgorithm::TTPCTRExPatSubAlgorithm(trex::TTPCTRExPatSubAlgorithm&& in){
  fHitMap=std::move(in.fHitMap);
  fVolGroupMan=in.fVolGroupMan;
  fAStar=in.fAStar;
  fTracks=std::move(in.fTracks);
  fHasHits=in.fHasHits;
  fHasValidPaths=in.fHasValidPaths;
  fTPC=in.fTPC;
  fLayout=in.fLayout;
  fPrimary=in.fPrimary;

  in.fAStar=0;
  in.fVolGroupMan=0;
}

void trex::TTPCTRExPatSubAlgorithm::SetUpHits(std::map<long, trex::TTPCUnitVolume*>& map, trex::TTPCAStar* aStarCopy){

  fHitMap=map;
  
  fHasHits = fHitMap.size() > 0;
  // add map of hits to group manager object
  fVolGroupMan->AddPrimaryHits(fHitMap);
  // add map of hits to path finder
  fAStar->AddHits(fVolGroupMan, fHitMap);
}

void trex::TTPCTRExPatSubAlgorithm::ProduceContainers(){

  std::cout<<"1"<<std::endl;
  // ignore if hits don't already exist
  if(!fHasHits) return;

  // start out having no valid paths by default
  fHasValidPaths = false;

  // find list of hit groups on edge of paths
  std::vector< trex::TTPCVolGroup > edgeGroups=fVolGroupMan->GetEdgeGroups();

  std::cout<<"Subalg finds "<<edgeGroups.size()<<" initial edge groups"<<std::endl;

  // clear redundant edge groups (i.e. ones which are actually just half way points along existing paths)
  fAStar->ClearRedundancies(fVolGroupMan, edgeGroups);

  std::cout<<"Subalg finds "<<edgeGroups.size()<<" non-redundant edge groups"<<std::endl;

  // connect pairs of edge groups to find missing hits
  std::vector< trex::TTPCOrderedVolGroup > edgePaths;
  fAStar->ConnectGroupsOrdered(fVolGroupMan, edgeGroups, edgePaths, false, false);

  // find unmatched hits
  std::map<long, trex::TTPCUnitVolume*> unmatchedVolumes;
  unmatchedVolumes = fHitMap;
  
  for(std::vector< trex::TTPCOrderedVolGroup>::iterator edgePathIt = edgePaths.begin(); edgePathIt != edgePaths.end(); ++edgePathIt){
    trex::TTPCOrderedVolGroup& edgePath = *edgePathIt;

    fVolGroupMan->BuildGroupFriends(edgePath, trex::TTPCConnection::extraHits);
    trex::TTPCVolGroup pathHits=edgePath.GetExtendedHits();

    for(std::map<long, trex::TTPCUnitVolume*>::iterator foundVolIt = pathHits.begin(); foundVolIt != pathHits.end(); ++foundVolIt){
      long id = foundVolIt->first;
      std::map<long, trex::TTPCUnitVolume*>::iterator foundIdLoc = unmatchedVolumes.find(id);
      if(foundIdLoc != unmatchedVolumes.end()) unmatchedVolumes.erase(foundIdLoc);
    };
  };

  // split and group unmatched hits
  trex::TTPCVolGroupMan secondPassVolume(fLayout);
  secondPassVolume.AddPrimaryHits(unmatchedVolumes);
  std::vector< trex::TTPCVolGroup > extraSubVolumes;
  secondPassVolume.GetConnectedHits(extraSubVolumes,trex::TTPCConnection::path, trex::TTPCHitGroupings::all, true);

  for(std::vector< trex::TTPCVolGroup >::iterator extraSubVolumeIt = extraSubVolumes.begin(); extraSubVolumeIt != extraSubVolumes.end(); ++extraSubVolumeIt){
    trex::TTPCVolGroup& extraSubVolume = *extraSubVolumeIt;

    trex::TTPCVolGroupMan subVolumeMan(fLayout);

    //MDH
    //Hmm, bit clumsy. Add a getter for the TTPCVolGroup hit map?
    std::map<long, trex::TTPCUnitVolume*> subVolumeHits;
    for(std::map<long, trex::TTPCUnitVolume*>::iterator volIt = extraSubVolume.begin(); volIt != extraSubVolume.end(); ++volIt){
      subVolumeHits[volIt->first] = volIt->second;
    }
    subVolumeMan.AddPrimaryHits(subVolumeHits);
    std::vector<trex::TTPCVolGroup> extraEdgeGroups=subVolumeMan.GetEdgeGroups();
    // add unmatched edges to main group
    for(std::vector< trex::TTPCVolGroup >::iterator extraEdgeGroupIt = extraEdgeGroups.begin(); extraEdgeGroupIt != extraEdgeGroups.end(); ++extraEdgeGroupIt){
      edgeGroups.emplace_back(std::move(*extraEdgeGroupIt));
    }
  }
  std::cout<<"Subalg finds "<<edgeGroups.size()<<" edge groups (2nd pass)"<<std::endl;
  // clear redundant edge groups (i.e. ones which are actually just half way points along existing paths)
  fAStar->ClearRedundancies(fVolGroupMan, edgeGroups);

  std::cout<<"Subalg finds "<<edgeGroups.size()<<" non-redundant edge groups (2nd pass)"<<std::endl;
  
  // connect pairs of edge groups

  // define container for true paths before corner detection
  std::vector< trex::TTPCOrderedVolGroup > truePaths;

  // if there are at least three edge groups, look for vertices to connect to them
  if(edgeGroups.size() > 2){

    edgePaths.clear();

    // connect pairs of edge groups
    fAStar->ConnectGroupsOrdered(fVolGroupMan, edgeGroups, edgePaths, false, true);
    // look for vertices from edge group pairs
    fVolGroupMan->GetFoci(edgePaths, fVertices ,.1);
    // make sure number of vertices is two less than number of track ends
    fVolGroupMan->CleanUpVertices(edgeGroups, fVertices);
    // make sure vertices contain all the right hits

    fVolGroupMan->BulkGroups(fVertices);

    if(fVertices.size() > 0){
      // connect vertices to path ends to get paths
      fAStar->ConnectVertexGroupsOrdered(fVolGroupMan, fVertices, edgeGroups,truePaths);
      // clear redundant paths
      fAStar->ClearVertexConnectionRedundancies(fVolGroupMan, truePaths, fVertices);
    }
    // otherwise, return empty handed
    else{
      return;
    };
  }
  // if there are two, connect the two edge groups
  else if (edgeGroups.size() > 1){
    //fAStar->ConnectGroupsOrdered(fVolGroupMan, edgeGroups, truePaths);
    truePaths=std::move(edgePaths);
  }
  // otherwise, return empty handed
  else{
    return;
  };

  std::cout<<"After joining to vertices have "<<truePaths.size()<<" groups"<<std::endl;
  
  // ensure that paths don't share hits with each other
  fVolGroupMan->BuildAllFriends(truePaths);

  // extend and cluster paths
  for(std::vector< trex::TTPCOrderedVolGroup >::iterator truePathIt = truePaths.begin(); truePathIt != truePaths.end(); ++truePathIt){
    fVolGroupMan->ClusterGroupFriends(*truePathIt, false, true);
  };
  // clear empties
  
  std::vector<trex::TTPCOrderedVolGroup> nonEmptyPaths;
  for(std::vector< trex::TTPCOrderedVolGroup >::iterator pathKeep = truePaths.begin();pathKeep != truePaths.end();++pathKeep){
  if(pathKeep->size()) nonEmptyPaths.emplace_back(std::move(*pathKeep));
  }

  truePaths=std::move(nonEmptyPaths);

  std::cout<<"After clearing empties have "<<truePaths.size()<<" groups"<<std::endl;

  // look for kinks in paths
  fVolGroupMan->BreakPathsAboutKinks(truePaths);

  for(std::vector< trex::TTPCOrderedVolGroup >::iterator brokenPathIt = truePaths.begin(); brokenPathIt != truePaths.end(); ++brokenPathIt){
    fVolGroupMan->ClusterGroupFriends(*brokenPathIt, true, true);
  };

  std::cout<<"After splitting kinks have "<<truePaths.size()<<" groups"<<std::endl;

  // temporarily merge x-paths into junctions to make the next two steps easier
  fVolGroupMan->SeparateXPathHits(truePaths);

  // clean up any and all instances of hits beiing placed in the wrong object
  fVolGroupMan->SeparateEmptyClusters(truePaths);
  fVolGroupMan->SeparateClusterHits(truePaths);
  fVolGroupMan->SeparateAnomHits(truePaths);
  fVolGroupMan->SeparateJunctionHits(truePaths);

  // enforce path ordering
  fVolGroupMan->EnforceOrdering(truePaths);
  std::cout<<"After splitting off hits have "<<truePaths.size()<<" groups"<<std::endl;
  // add all paths that make sense and store if none are found
  fVolGroupMan->SanityFilter(truePaths);
  std::cout<<"After sanity filter have "<<truePaths.size()<<" groups"<<std::endl;
  // associate any remaining unused hits to junctions where possible
  trex::TTPCVolGroup unusedHits(fLayout);
  
  std::cout << "DEBUGGING: Accessing unused Hits" << std::endl;

  fVolGroupMan->GetUnusedHits(truePaths,unusedHits);
  fVolGroupMan->AssociateUnusedWithJunctions(unusedHits, truePaths);

  std::cout << "DEBUGGING: Resetting vertex status" << std::endl;

  // reset vertex status based on how many paths a junction is connected to
  fVolGroupMan->ResetVertexStatuses(truePaths);

  std::cout << "DEBUGGING: moving truePaths to fTracks" << std::endl;
  

  fTracks=std::move(truePaths);

  std::cout << "DEBUGGING: checking has valid paths" << std::endl;

  fHasValidPaths = fTracks.size() > 0;


  std::cout << "DEBUGGING: Produce Containers has run!" << std::endl;
}

//trex::TTPCPattern* trex::TTPCTRExPatSubAlgorithm::GetPattern(){
//  if(!fPattern){
//    return 0;
//  };
//  return fPattern;
//}

void trex::TTPCTRExPatSubAlgorithm::ProducePattern(TTRExPattern& output){//trex::THitSelection* used){

  std::cout << "DEBUGGING: 1" << std::endl;

  std::vector< std::vector<unsigned int> > junctionsToPathsMap;
  std::vector<trex::TTRExPath>& paths=output.GetPaths();
  std::vector<trex::TTRExJunction>& juncts=output.GetJunctions();


  std::cout << "DEBUGGING: 2" << std::endl;

  std::vector< trex::TTPCVolGroup* > junctionGroups;

  std::cout << "DEBUGGING: 3" << std::endl;

  for(std::vector< trex::TTPCOrderedVolGroup >::iterator trackIt = fTracks.begin(); trackIt != fTracks.end(); ++trackIt){
    trex::TTPCOrderedVolGroup& track = *trackIt;

    // ignore if empty
    std::vector<trex::TTRExHVCluster*> clusters = track.GetClusters();
    if(!clusters.size()) continue;

    paths.emplace_back(clusters);

    // and save if the path is in the x direction
    //    if(track->GetIsXPath()){
    //  path->SetIsXPath(true);
    //};

    // save groups for junctions
    std::vector< trex::TTPCVolGroup* > pathJunctionGroups;
    if(track.HasBackHits()) pathJunctionGroups.push_back(&(track.GetBackHits()));
    if(track.HasFrontHits()) pathJunctionGroups.push_back(&(track.GetFrontHits()));
    
    std::cout << "DEBUGGING: 4" << std::endl;

    // determine if junction candidate has already been found and create a group for it if it hasn't
    for(std::vector< trex::TTPCVolGroup* >::iterator junctionGroupIt = pathJunctionGroups.begin(); junctionGroupIt != pathJunctionGroups.end(); ++junctionGroupIt){
      trex::TTPCVolGroup* junctionGroup = *junctionGroupIt;
      bool found = false;

      int iMax = junctionGroups.size();
      for(int i=0; i<iMax; i++){

	std::cout << "DEBUGGING: 5" << std::endl;

        // if the junction has the same id as the temporary, add this path to it
        if(junctionGroup->GetID() == junctionGroups[i]->GetID()){
          // if already found, check this path isn't already added and add this path
          bool pathInJunction = false;
          for(auto constituentIt = junctionsToPathsMap[i].begin(); 
	      constituentIt != junctionsToPathsMap[i].end();
	      ++constituentIt){
            if(paths.size()-1 == *constituentIt){
              pathInJunction = true;
              break;
            };
          };
          if(!pathInJunction) junctionsToPathsMap[i].push_back(paths.size()-1);

	  std::cout << "DEBUGGING: 5" << std::endl;

          found = true;
        };
      };


      // if not alread found, create a new one
      if(!found){


	/*        TVector3 tempVect3 = junctionGroup->GetAverageVol()->GetPos();
        double tempTime = junctionGroup->GetAverageVol()->GetTime();
        // Get a guess X position based on a default T0 of 0.0 just to have a useful X position downstream.
        double tempX = TTPCUtils::GetXWithT0(tempTime, 0.0, (*(junctionHitSelection->begin())));
        tempVect3.SetX(tempX);
        TLorentzVector tempPosition(tempVect3, tempTime);
	*/


	//	juncts.emplace_back(junctionGroup->GetHits());
        junctionsToPathsMap.emplace_back();
	junctionsToPathsMap.back().push_back(paths.size()-1);
	junctionGroups.push_back(junctionGroup);

      };
    };
  };

  std::cout << "DEBUGGING: 6" << std::endl;

  // Connect Paths and Junctions together according to the Map
  for(int i=0; i<junctionsToPathsMap.size(); ++i){
    std::cout << "bug 1" << std::endl;
    

    std::cout << "Size of Map " << junctionsToPathsMap.size() << " Compared to Size of JunctionGoups " << junctionGroups.size() << std::endl;


    //WHY?
    //if(junctionsToPathsMap[i].size()<2) continue;
    
    std::cout << "bug 2" << std::endl;
    
    juncts.emplace_back(junctionGroups[i]->GetHits());

    std::cout << "bug 3" << std::endl;

    for(int j=0; j<junctionsToPathsMap[i].size(); ++j){

      std::cout << "bug 4" << std::endl;
      int pathIndex = junctionsToPathsMap[i][j];
      std::cout << "bug 5" << std::endl;

      std::cout << "Length of path: " << paths.size() << " Compared to path Index " << pathIndex << std::endl;
      std::cout << "Length of juncts: " << juncts.size() << " Compared to junct Index " << i << std::endl;

      ConnectJunctionAndPath(juncts[i], paths[pathIndex]);
      std::cout << "bug 6" << std::endl;
      
    }
    std::cout << "bug 7" << std::endl;
  }

  std::cout << "DEBUGGING: 7" << std::endl;

  std::cout<<"Built a pattern..."<<std::endl;
  std::cout<<"  from "<<fHitMap.size()<<" available hits"<<std::endl;
  std::cout<<"  "<<paths.size()<<" paths"<<std::endl;
  std::cout<<"  and  "<<juncts.size()<<" junctions"<<std::endl;
  for(int i=0;i<paths.size();++i){
    std::vector<trex::TTRExHVCluster*> clusters = paths[i].GetClusters();

    int nHits=0;
    for(auto iCl=clusters.begin();iCl!=clusters.end();++iCl){nHits+=(*iCl)->GetClusterHits().size();}

    std::cout<<"   Path "<<i<<" has "<<clusters.size()<<" clusters"<<std::endl;
    std::cout<<"   and "<<nHits<<" hits"<<std::endl;
    std::cout<<"  **********"<<std::endl;
  }
  for(int i=0;i<juncts.size();++i){
    std::cout<<"   Junction "<<i<<" has "<<juncts[i].GetHits().size()<<" hits"<<std::endl;
    std::cout<<"   and is linked to paths:";
    for(int j=0;j<junctionsToPathsMap[i].size();++j){
      std::cout<<(j==0?" ":", ")<<junctionsToPathsMap[i][j]<<std::endl;
    }
    std::cout<<"  **********"<<std::endl;
  }
  std::cout<<"**************"<<std::endl;
  // fill used and unused hits
  //  FillUsedHits(used);

}

std::vector<trex::TTPCHitPad*> trex::TTPCTRExPatSubAlgorithm::GetHits(){
  return fVolGroupMan->GetHits();
}


std::vector<trex::TTRExHVCluster*> trex::TTPCTRExPatSubAlgorithm::GetHits(trex::TTPCOrderedVolGroup& path){
  return std::move(path.GetClusters());
}

void trex::TTPCTRExPatSubAlgorithm::GetRegions(std::vector< trex::TTPCVolGroup >& regions){

  if(!fHasHits){
    regions.clear();
    return;
  }
  // get list of hit groups for sub-events
  fVolGroupMan->GetConnectedHits(regions,trex::TTPCConnection::path);
}

//MDH
//Not currently used. May be used later if we resurrect ProducePattern.
/*
void trex::TTPCTRExPatSubAlgorithm::FillUsedHits(std::vector<trex::TTPCHitPad*>& used){
  // if pattern hasn't been correctly built, hits can't possibly be used
  if(!fPattern) return;

  // matched (initially empty) volumes
  std::map<long, trex::TTPCUnitVolume*> matchedVolumes = std::map<long, trex::TTPCUnitVolume*>();

  // find those contained in any paths or junctions
  for(std::vector< trex::THandle<trex::TTPCOrderedVolGroup> >::iterator trackIt = fTracks.begin(); trackIt != fTracks.end(); ++trackIt){
    trex::THandle<trex::TTPCOrderedVolGroup> track = *trackIt;
    // path hits
    for(std::vector<trex::TTPCPathVolume*>::iterator pathIt = track->begin(); pathIt != track->end(); ++pathIt){
      trex::TTPCPathVolume* path = *pathIt;

      std::vector<trex::TTPCUnitVolume*> pathClusterVols = path->GetExtendedCell(0);
      for(std::vector<trex::TTPCUnitVolume*>::iterator volIt = pathClusterVols.begin(); volIt != pathClusterVols.end(); ++volIt){
        std::map<long, trex::TTPCUnitVolume*>::iterator foundVol = fHitMap.find((*volIt)->GetID());
        matchedVolumes[foundVol->first] = foundVol->second;
      };
    };
  };

  // junction hits
  std::vector< trex::THandle<trex::TTPCVolGroup> > pathVertices = fVolGroupMan->GetJunctionsFromPaths(fTracks);
  for(std::vector< trex::THandle<trex::TTPCVolGroup> >::iterator pathVertexIt = pathVertices.begin(); pathVertexIt != pathVertices.end(); ++pathVertexIt){
    trex::THandle<trex::TTPCVolGroup> pathVertex = *pathVertexIt;

    for(std::map<long, trex::TTPCUnitVolume*>::iterator volIt = pathVertex->begin(); volIt != pathVertex->end(); ++volIt){
      std::map<long, trex::TTPCUnitVolume*>::iterator foundVol = fHitMap.find(volIt->first);
      matchedVolumes[foundVol->first] = foundVol->second;
    };
  };

  // add used hits
  for(std::map<long, trex::TTPCUnitVolume*>::iterator volEl = matchedVolumes.begin(); volEl != matchedVolumes.end(); ++volEl){
    trex::TTPCUnitVolume* vol = volEl->second;
    if(!vol) continue;

    std::vector< trex::THandle<trex::TTPCHitPad> > hits = vol->GetHits();
    for(std::vector< trex::THandle<trex::TTPCHitPad> >::iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt){
      trex::THandle<trex::TTPCHitPad> hitPad = *hitIt;
      if(hitPad) used->AddHit(hitPad);
    };
  };
}
*/
void trex::TTPCTRExPatSubAlgorithm::GetTrackExtendedHits(std::vector<trex::TTPCVolGroup>& extHits){

  for(std::vector< trex::TTPCOrderedVolGroup >::iterator trackIt = fTracks.begin(); trackIt != fTracks.end(); ++trackIt){
    trex::TTPCOrderedVolGroup& track = *trackIt;
    extHits.push_back(track.GetExtendedHits());
  }
  
}


void trex::TTPCTRExPatSubAlgorithm::ConnectJunctionAndPath(trex::TTRExJunction& junction, trex::TTRExPath& path){

  trex::TTRExPath *pat = &path;
  trex::TTRExJunction *junct = &junction;

  junction.AddConnectedPath(pat);
  path.AddConnectedJunction(junct);
}
