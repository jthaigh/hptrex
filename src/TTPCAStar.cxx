// eddy
#include "TTPCAStar.hxx"

trex::TTPCAStar::TTPCAStar(trex::TTPCLayout* layout){
  fLayout = layout;

  fAStarPoints = std::vector<trex::TTPCAStarPoint*>();
  fHitMap = std::map<long, trex::TTPCUnitVolume*>();

  // defualt scales and heuristic factor
  fXScale = fLayout->GetAStarXScale();
  fYScale = fLayout->GetAStarYScale();
  fZScale = fLayout->GetAStarZScale();

  fHeuristicFactor = fLayout->GetAStarHeuristicFactor();
}
trex::TTPCAStar::~TTPCAStar(){
  // delete all A* points
  for(std::vector<trex::TTPCAStarPoint*>::iterator pnt = fAStarPoints.begin(); pnt != fAStarPoints.end(); ++pnt){
    delete *pnt;
  };
}

void trex::TTPCAStar::AddHits(trex::TTPCVolGroupMan* volGroupMan, std::map<long, trex::TTPCUnitVolume*> hitMap, bool extendMode){
  // add hit map
  fHitMap = hitMap;

  // set up vector of points
  fAStarPoints = std::vector<trex::TTPCAStarPoint*>();

  // add all hits to list
  for(std::map<long, trex::TTPCUnitVolume*>::iterator hitEl = hitMap.begin(); hitEl != hitMap.end(); hitEl++){
    trex::TTPCAStarPoint* pnt = new trex::TTPCAStarPoint();

    pnt->vol = hitEl->second;

    pnt->x = hitEl->second->GetX();
    pnt->y = hitEl->second->GetY();
    pnt->z = hitEl->second->GetZ();

    pnt->aStarFriends = std::map<trex::TTPCAStarPoint*, float>();
    pnt->aStarHasFriends = false;

    fAStarPoints.push_back(pnt);
  };

  for(std::vector<trex::TTPCAStarPoint*>::iterator pointIt = fAStarPoints.begin(); pointIt != fAStarPoints.end(); ++pointIt){
    trex::TTPCAStarPoint* point = *pointIt;
    // build connections now and share friends backwards to avoid asymmetries
    GetNearHitConnections(volGroupMan, point, extendMode);
    for(std::map<trex::TTPCAStarPoint*, float>::iterator fr = point->aStarFriends.begin(); fr != point->aStarFriends.end(); ++fr){
      fr->first->aStarFriends[point] = fr->second;
    };
  };
}

void trex::TTPCAStar::AddHits(trex::TTPCAStar* prevAStar, std::map<long, trex::TTPCUnitVolume*> hitMap){
  // add hit map
  fHitMap = hitMap;

  // set up vector of points
  fAStarPoints = std::vector<trex::TTPCAStarPoint*>();

  // map to help copying points over
  std::map<trex::TTPCAStarPoint*, trex::TTPCAStarPoint*> oldNewMap;
  // copy hits that are in both the previous A* container and the map
  for(std::vector<trex::TTPCAStarPoint*>::iterator oldPntIt = prevAStar->begin(); oldPntIt != prevAStar->end(); ++oldPntIt){
    // make sure point is inside map
    trex::TTPCAStarPoint* oldPnt = *oldPntIt;
    long id = oldPnt->vol->GetID();

    if(hitMap.find(id) != hitMap.end()){
      trex::TTPCAStarPoint* pnt = new trex::TTPCAStarPoint();

      pnt->vol = oldPnt->vol;

      pnt->x = oldPnt->x;
      pnt->y = oldPnt->y;
      pnt->z = oldPnt->z;

      pnt->aStarFriends = std::map<trex::TTPCAStarPoint*, float>();

      oldNewMap[oldPnt] = pnt;
    }
  };

  // now add friends
  for(std::vector<trex::TTPCAStarPoint*>::iterator oldPntIt = prevAStar->begin(); oldPntIt != prevAStar->end(); ++oldPntIt){
    trex::TTPCAStarPoint* oldPnt = *oldPntIt;
    long id = oldPnt->vol->GetID();

    if(hitMap.find(id) != hitMap.end()){
      trex::TTPCAStarPoint* pnt = oldNewMap.at(oldPnt);
      // copy valid friends over too
      for(std::map<trex::TTPCAStarPoint*, float>::iterator fr = oldPnt->aStarFriends.begin(); fr != oldPnt->aStarFriends.end(); ++fr){
        long id = oldPnt->vol->GetID();
        if(hitMap.find(id) != hitMap.end()){
          pnt->aStarFriends[ oldNewMap.at(fr->first) ] = fr->second;
        };
      };

      pnt->aStarHasFriends = true;

      fAStarPoints.push_back(pnt);
    };
  };
}

//MDH
//Not used
/*
std::vector< trex::THandle<trex::TTPCVolGroup> > trex::TTPCAStar::ConnectVertexGroups(trex::TTPCVolGroupMan* volGroupMan, std::vector< trex::THandle<trex::TTPCVolGroup> > vertices, std::vector< trex::THandle<trex::TTPCVolGroup> > edges, int maxNo){
  std::vector< trex::THandle<trex::TTPCVolGroup> > connections;
  // find ordered sets of connections between the groups
  std::vector< trex::THandle<trex::TTPCOrderedVolGroup> > orderedConnections = ConnectVertexGroupsOrdered(volGroupMan, vertices, edges, maxNo);

  // convert the ordered groups to unordered groups
  for(std::vector< trex::THandle<trex::TTPCOrderedVolGroup> >::iterator path = orderedConnections.begin(); path != orderedConnections.end(); ++path){
    trex::THandle<trex::TTPCVolGroup> connection(new trex::TTPCVolGroup(fLayout));
    for(std::vector<trex::TTPCPathVolume*>::iterator cell = (*path)->begin(); cell != (*path)->end(); ++cell){
      connection->AddHit((*cell)->GetUnitVolume());
    };
    connections.push_back(connection);
  };

  return connections;
}
*/

void trex::TTPCAStar::ConnectVertexGroupsOrdered(trex::TTPCVolGroupMan* volGroupMan, std::vector< trex::TTPCVolGroup >& vertices, std::vector< trex::TTPCVolGroup >& edges, std::vector< trex::TTPCOrderedVolGroup >& connections, int maxNo){

  // keep track of number of iterations to break when exceeded
  int i=0;
  // loop over each vertex and each edge
  for(std::vector< trex::TTPCVolGroup >::iterator ver = vertices.begin(); ver != vertices.end(); ++ver){
    for(std::vector< trex::TTPCVolGroup >::iterator grp = edges.begin(); grp != edges.end(); ++grp){
      // add connection between the groups to set of connections
      connections.emplace_back(fLayout);
      ConnectGroupPair(*ver,*grp,connections.back(), true,false);
      if(connections.back().empty()) connections.pop_back();
      // break out if the number of iterations is too high
      i++;
      if(i>=maxNo) break;
    };
    if(i>=maxNo) break;
  };

  ConnectGroupsOrdered(volGroupMan, vertices, connections, true, false, maxNo);

}

//MDH
//Not used
/*
std::vector< trex::THandle<trex::TTPCVolGroup> > trex::TTPCAStar::ConnectGroups(trex::TTPCVolGroupMan* volGroupMan, std::vector< trex::THandle<trex::TTPCVolGroup> > groups, bool allConnections, int maxNo){
  std::vector< trex::THandle<trex::TTPCVolGroup> > connections;
  // find ordered sets of connections between the groups
  std::vector< trex::THandle<trex::TTPCOrderedVolGroup> > orderedConnections = ConnectGroupsOrdered(volGroupMan, groups, false, allConnections, maxNo);

  // convert the ordered groups to unordered groups
  for(std::vector< trex::THandle<trex::TTPCOrderedVolGroup> >::iterator path = orderedConnections.begin(); path != orderedConnections.end(); ++path){
    trex::THandle<trex::TTPCVolGroup> connection(new trex::TTPCVolGroup(fLayout));
    for(std::vector<trex::TTPCPathVolume*>::iterator cell = (*path)->begin(); cell != (*path)->end(); ++cell){
      connection->AddHit((*cell)->GetUnitVolume());
    };
    connections.push_back(connection);
  };

  return connections;
}
*/

void trex::TTPCAStar::ConnectGroupsOrdered(trex::TTPCVolGroupMan* volGroupMan, std::vector< trex::TTPCVolGroup >& groups, std::vector< trex::TTPCOrderedVolGroup >& connections, bool vertices, bool allConnections, int maxNo){

  // keep track of number of iterations to break when exceeded
  int i=0;
  // loop over each pair of groups once
  for(std::vector< trex::TTPCVolGroup >::iterator grp1 = groups.begin(); grp1 != groups.end(); ++grp1){
    std::vector< trex::TTPCVolGroup >::iterator grp2Begin = allConnections ? groups.begin() : grp1+1;
    std::vector< trex::TTPCVolGroup >::iterator grp2End = groups.end();
    for(std::vector< trex::TTPCVolGroup >::iterator grp2 = grp2Begin; grp2 != grp2End; ++grp2){
      if(grp1 == grp2) continue;
      i++;
      // add connection between the groups to set of connections
      connections.emplace_back(fLayout);
      ConnectGroupPair(*grp1,*grp2, connections.back(), vertices,vertices);
      if(connections.back().empty()) connections.pop_back();
      // break out if the number of iterations is too high
      if(i>=maxNo) break;
    };
    if(i>=maxNo) break;
  };

}

//Ok this needs a bit of work. For now only changed the signature
//so that this will be an in-place cleanup
void trex::TTPCAStar::ClearRedundancies(trex::TTPCVolGroupMan* volGroupMan, std::vector< trex::TTPCVolGroup >& groups, int maxNo){
  // the next block of code is pretty ugly but needed to remove bugs in some very specific topologies
  // first, merge groups which are too close together

  std::vector< trex::TTPCVolGroup > mergedGroups;

  std::cout<<"ClearRedundancies starting with "<<groups.size()<<" groups"<<std::endl;

  for(auto iGrp=groups.begin();iGrp!=groups.end();++iGrp){
    TVector3 pos=iGrp->GetAveragePosition();
    std::cout<<"Group at "<<pos.X()<<", "<<pos.Y()<<", "<<pos.Z()<<std::endl;
  }

  while(groups.size()){
    // copy id of merge from first inGroup for consistency
    unsigned int id = (*groups.begin()).GetID();

    std::vector< std::vector< trex::TTPCVolGroup >::iterator > mergeIts;
    mergeIts.push_back(groups.begin());

    // stop only when all groups overlapping are merged
    bool furtherMerges = true;
    while(furtherMerges){
      furtherMerges = false;

      // add all groups overlapping with one in the 'to merge' list
      if(groups.size() > 1)
      for(std::vector< std::vector< trex::TTPCVolGroup >::iterator >::iterator mergeItIt = mergeIts.begin(); mergeItIt != mergeIts.end(); ++mergeItIt){
        std::vector< trex::TTPCVolGroup >::iterator mergeIt = *mergeItIt;
        trex::TTPCVolGroup& merge = *mergeIt;

        for(std::vector< trex::TTPCVolGroup >::iterator inGroupIt = groups.begin(); inGroupIt != groups.end(); ++inGroupIt){
          trex::TTPCVolGroup& inGroup = *inGroupIt;

          bool alreadyAdded = false;
          for(std::vector< std::vector< trex::TTPCVolGroup >::iterator >::iterator mergeItIt2 = mergeIts.begin(); mergeItIt2 != mergeIts.end(); ++mergeItIt2)
          if(*mergeItIt2 == inGroupIt){
            alreadyAdded = true;
            break;
          };
          if(alreadyAdded) continue;

          if( volGroupMan->GetGroupGroupOverlap(inGroup, merge, trex::TTPCConnection::edgeMerge, true) ){
            mergeIts.push_back(inGroupIt);
            furtherMerges = true;
            if(furtherMerges) break;
          };
        };
        if(furtherMerges) break;
      };
    };
  
    // set new group with the ID of the first of the merges
    mergedGroups.emplace_back(fLayout,id);
    trex::TTPCVolGroup& groupOut = mergedGroups.back();
    std::vector< std::vector<trex::TTPCVolGroup>::iterator > delGroups;

    for(std::vector< std::vector< trex::TTPCVolGroup >::iterator >::iterator mergeItIt = mergeIts.begin(); mergeItIt != mergeIts.end(); ++mergeItIt){
      std::vector< trex::TTPCVolGroup >::iterator mergeIt = *mergeItIt;
      trex::TTPCVolGroup& merge = *mergeIt;

      groupOut.AddHits(merge);
      groupOut.SetXLean(merge.GetXLean());
      groupOut.SetYLean(merge.GetYLean());
      groupOut.SetZLean(merge.GetZLean());

      delGroups.push_back(mergeIt);
    };

    std::vector< trex::TTPCVolGroup > newInGroups;
    for(std::vector< trex::TTPCVolGroup >::iterator grpIt = groups.begin(); grpIt != groups.end(); ++grpIt){
      if(std::find(delGroups.begin(),delGroups.end(),grpIt)==delGroups.end()){
      newInGroups.push_back(*grpIt);
      }
    }
      groups = std::move(newInGroups);
  }

  std::cout<<"After removing overlaps have "<<mergedGroups.size()<<" groups"<<std::endl;

  // parallel arrays for groups and list of groups they overlap
  std::vector< std::vector<int> > connectedGroups;
  for(std::vector< trex::TTPCVolGroup >::iterator grpIt = mergedGroups.begin(); grpIt != mergedGroups.end(); ++grpIt){
    connectedGroups.emplace_back();
  }
  std::vector<bool> groupsMaySurvive (mergedGroups.size(), true);

  // keep track of number of iterations to break when exceeded
  int nGroups = (int)mergedGroups.size();
  int count=0;

  // loop over each pair of groups once
  for(int i=0; i<nGroups; ++i){
    for(int j=i+1; j<nGroups; ++j){
      trex::TTPCVolGroup& grp1 = mergedGroups[i];
      trex::TTPCVolGroup& grp2 = mergedGroups[j];
      // add connection between the groups to set of connections
      trex::TTPCOrderedVolGroup path(fLayout);
      ConnectGroupPair(grp1, grp2,path);
      // kill any groups associated with this path other than its direct ends
      for(int k=0; k<nGroups; ++k){
        if(k == i) continue;
        if(k == j) continue;
        trex::TTPCVolGroup& grp3 = mergedGroups[k];
        if(volGroupMan->GetPathVolOverlap(path, grp3.GetAverageVol())){
          if(std::find(connectedGroups[i].begin(), connectedGroups[i].end(), k) == connectedGroups[i].end()) connectedGroups[i].push_back(k);
          if(std::find(connectedGroups[j].begin(), connectedGroups[j].end(), k) == connectedGroups[j].end()) connectedGroups[j].push_back(k);
        };
      };
      // break out if the number of iterations is too high
      count++;
      if(count>=maxNo) break;
    };
    if(count>=maxNo) break;
  };


  // build surviving groups
  std::vector< trex::TTPCVolGroup > survivors;

  // first, add all totally safe groups
  for(int i=0; i<nGroups; ++i){
    bool madeRedundant = false;

    // is this group marked for death?
    if(!groupsMaySurvive[i]){
      madeRedundant = true;
    }
    else{
      // is this group made redundant by any other?
      for(int j=0; j<nGroups; ++j){
        if(j==i) continue;

        int connectionSize = (int)connectedGroups[j].size();
        for(int k=0; k<connectionSize; ++k){
          if(i==connectedGroups[j][k]) madeRedundant = true;
          if(madeRedundant) break;
        if(madeRedundant) break;
        };
      };
    };

    if(!madeRedundant){
      survivors.push_back(mergedGroups[i]);

      // mark everything this group overlapped with for death
      int connection2Size = connectedGroups[i].size();
      for(int conn=0; conn<connection2Size; ++conn) groupsMaySurvive[ connectedGroups[i][conn] ] = false;
      connectedGroups[i].clear();
    }
  }

  // look for loops and merge them
  for(int i=0; i<nGroups; ++i){
    // only include if isolated
    if(!groupsMaySurvive[i]) continue;

    int connectionSize = connectedGroups[i].size();
    for(int j=0; j<connectionSize; j++){
      int merge1 = i;
      int merge2 = connectedGroups[i][j];
      // only include if isolated
      if(!groupsMaySurvive[merge2]) continue;

      // avoid duplication
      if(merge2 > merge1){
        int connectionMergeSize = (int)connectedGroups[merge2].size();
        for(int k=0; k<connectionMergeSize; k++){
          if(connectedGroups[merge2][k] == merge1){
            // first, check if either is made redundant
	    survivors.emplace_back(fLayout,mergedGroups[merge1].GetID());
            volGroupMan->MergeGroups(mergedGroups[merge1], mergedGroups[merge2], survivors.back());
          };
          break;
        };
      };
    };
  };
  std::cout<<"After removing midpoint groups have "<<survivors.size()<<" groups"<<std::endl;
  groups=std::move(survivors);
}

void trex::TTPCAStar::ClearVertexConnectionRedundancies(trex::TTPCVolGroupMan* volGroupMan, std::vector< trex::TTPCOrderedVolGroup >& paths, std::vector< trex::TTPCVolGroup >& vertices){
  std::vector< trex::TTPCOrderedVolGroup > outPaths;

  for(std::vector< trex::TTPCOrderedVolGroup >::iterator path = paths.begin(); path != paths.end(); ++path){
    bool redundancy = false;
    for(std::vector< trex::TTPCVolGroup >::iterator vert = vertices.begin(); vert != vertices.end(); ++vert){
      if(volGroupMan->GetPathVolOverlap(*path, vert->GetAverageVol(), trex::TTPCConnection::vertexPath)){
	redundancy=true;
        break;
      }
    }
    if(!redundancy){
      outPaths.emplace_back(std::move(*path));
    }
  }

  paths=std::move(outPaths);
  
}

void trex::TTPCAStar::ConnectGroupPair(trex::TTPCVolGroup& group1, trex::TTPCVolGroup& group2, trex::TTPCOrderedVolGroup& connection, bool vertexGroup1, bool vertexGroup2, bool extendMode){

  // get unique ids and A* indices
  trex::TTPCUnitVolume* startCell = group1.GetAverageVol();
  trex::TTPCUnitVolume* endCell = group2.GetAverageVol();

  if(!startCell || !endCell) return;

  // reset relevant variables each time a connection needs to be made
  RebootHits();

  // assign A* points based on these
  trex::TTPCAStarPoint* startPoint=0;
  trex::TTPCAStarPoint* endPoint=0;
  for(std::vector<trex::TTPCAStarPoint*>::iterator pnt = fAStarPoints.begin(); pnt != fAStarPoints.end(); ++pnt){
    if((*pnt)->vol == startCell) startPoint = *pnt;
    if((*pnt)->vol == endCell) endPoint = *pnt;
  };

  // set up the A* variables chaining start point to end point
  int connectionReturnCode = DoConnection(startPoint, endPoint, extendMode);
  if(connectionReturnCode){
    return;
  };

  // add front and back groups to the returned group (front corresponds to hits near size()-1 index, back to hits near 0 index)
  connection.AddFrontHits(group1);
  connection.AddBackHits(group2);
  // set status of front and back groups as vertices (or not)
  connection.SetFrontIsVertex(vertexGroup1);
  connection.SetBackIsVertex(vertexGroup2);

  // start at end index and add the chain of cells connecting it from the start to the connection
  trex::TTPCAStarPoint* curPoint = endPoint;
  int nPoints = fAStarPoints.size();
  for(int i=0; i < nPoints; i++){
    connection.AddCell(curPoint->vol);

    trex::TTPCAStarPoint* nextPoint = curPoint->aStarParent;
    if(!nextPoint) break; 
    curPoint = nextPoint;
  };
  // add the first cell so long as the connection isn't one cell long
  if(curPoint != endPoint) connection.AddCell(curPoint->vol);

}

float trex::TTPCAStar::FindConnectionCost(trex::TTPCVolGroup& group1, trex::TTPCVolGroup& group2, bool extendMode, bool reduced, float maxCost){
  trex::TTPCUnitVolume* vol1 = group1.GetAverageVol();
  trex::TTPCUnitVolume* vol2 = group2.GetAverageVol();

  return FindConnectionCost(vol1, vol2, extendMode, reduced, maxCost);
}
float trex::TTPCAStar::FindConnectionCost(trex::TTPCVolGroup& group, trex::TTPCUnitVolume* vol, bool extendMode, bool reduced, float maxCost){
  trex::TTPCUnitVolume* groupVol = group.GetAverageVol();

  return FindConnectionCost(groupVol, vol, extendMode, reduced, maxCost);
}
float trex::TTPCAStar::FindConnectionCost(trex::TTPCUnitVolume* vol1, trex::TTPCUnitVolume* vol2, bool extendMode, bool reduced, float maxCost){
  float cost = -1.;
  if(!vol1 || !vol2) return cost;

  // reset relevant variables each time a connection needs to be made
  RebootHits();

  // assign A* points based on these
  trex::TTPCAStarPoint* startPoint=0;
  trex::TTPCAStarPoint* endPoint=0;
  for(std::vector<trex::TTPCAStarPoint*>::iterator pnt = fAStarPoints.begin(); pnt != fAStarPoints.end(); ++pnt){
    if((*pnt)->vol == vol1) startPoint = *pnt;
    if((*pnt)->vol == vol2) endPoint = *pnt;
  };

  int connectionReturnCode = DoConnection(startPoint, endPoint, extendMode, maxCost);
  if(connectionReturnCode){
    return cost;
  }
  else{
    cost = endPoint->aStarCost;

    if(cost <= 0.){
      cost = maxCost;
    };
  };

  return cost;
}

void trex::TTPCAStar::AssociateBestHits(trex::TTPCVolGroupMan* volGroupMan, std::vector< trex::TTPCOrderedVolGroup >& inPaths, bool extendMode, float maxCost){
  // add hits in path to group, making sure the same isn't added twice
  std::vector< trex::TTPCVolGroup > hitsGroups;
  std::set<trex::TTPCUnitVolume*> volsAdded;
  for(std::vector< trex::TTPCOrderedVolGroup >::iterator inPathIt = inPaths.begin(); inPathIt != inPaths.end(); ++inPathIt){
    trex::TTPCOrderedVolGroup& inPath = *inPathIt;

    hitsGroups.emplace_back(fLayout);
    trex::TTPCVolGroup& hitsGroup = hitsGroups.back();
    for(std::vector<trex::TTPCPathVolume*>::iterator pathVolIt = inPath.begin(); pathVolIt != inPath.end(); ++pathVolIt){
      trex::TTPCPathVolume* pathVol = *pathVolIt;
      trex::TTPCUnitVolume* vol = pathVol->GetUnitVolume();

      if(!volsAdded.count(vol)){
        hitsGroup.AddHit(vol);
        volsAdded.insert(vol);
      };
    };

  };

  // merge best hits in groups
  MergeBestHits(volGroupMan, hitsGroups, extendMode, maxCost);

  // add merged hits to appropriate paths
  for(unsigned int i=0; i<inPaths.size(); ++i){
    trex::TTPCOrderedVolGroup& inPath = inPaths.at(i);
    trex::TTPCVolGroup& hitsGroup = hitsGroups.at(i);

    inPath.AddExtendedHits(hitsGroup);
  };
}

//MDH
//Not used
/*
void trex::TTPCAStar::MergeBestHits(trex::TTPCVolGroupMan* volGroupMan, trex::THandle<trex::TTPCVolGroup> inGroup, bool fullASICPenalty, bool extendMode, float maxCost){
  std::vector< trex::THandle<trex::TTPCVolGroup> > dummyVector;
  dummyVector.push_back(inGroup);

  MergeBestHits(volGroupMan, dummyVector, fullASICPenalty, extendMode, maxCost);
}
*/

void trex::TTPCAStar::MergeBestHits(trex::TTPCVolGroupMan* volGroupMan, std::vector< trex::TTPCVolGroup >& inGroups, bool extendMode, float maxCost){
  // reset hits
  RebootHits();

  // define container of all near hits currently under consideration, starting with path hits
  std::vector<trex::TTPCAStarPoint*> openSet;
  std::vector<trex::TTPCAStarPoint*> definingSet (inGroups.size(), 0);
  int nGroup = 0;
  for(std::vector< trex::TTPCVolGroup >::iterator groupIt = inGroups.begin(); groupIt != inGroups.end(); ++groupIt){
    trex::TTPCVolGroup& group = *groupIt;

    for(std::map<long, trex::TTPCUnitVolume*>::iterator hitEl = group.begin(); hitEl != group.end(); ++hitEl){
      trex::TTPCUnitVolume* vol = hitEl->second;

      for(std::vector<trex::TTPCAStarPoint*>::iterator pntIt = fAStarPoints.begin(); pntIt != fAStarPoints.end(); ++pntIt){
        trex::TTPCAStarPoint* pnt = *pntIt;

        if(pnt->vol == vol){
          if(!definingSet[nGroup]){
            definingSet[nGroup] = pnt;
          };
          pnt->aStarParent = definingSet[nGroup];
          pnt->aStarCost = 0.;
          pnt->aStarOpen = true;
          openSet.push_back(pnt);
        };
      };
    };

    nGroup ++;
  };

  // loop over increasing max cost until either maximum distance is reached or no hits are left to add
  for(float i=1.; i<maxCost; i+=1.){
    // break out if all possible hits are checked
    if(!openSet.size()) break;

    // copy current open set over from old one
    std::vector<trex::TTPCAStarPoint*> curOpenSet;
    for(std::vector<trex::TTPCAStarPoint*>::iterator pntIt = openSet.begin(); pntIt != openSet.end(); ++pntIt){
      trex::TTPCAStarPoint* pnt = *pntIt;

      if(!pnt->aStarClosed){
        curOpenSet.push_back(*pntIt);
      };
    };
    openSet.clear();

    // loop once over current set
    for(std::vector<trex::TTPCAStarPoint*>::iterator pntIt = curOpenSet.begin(); pntIt != curOpenSet.end(); ++pntIt){
      trex::TTPCAStarPoint* pnt = *pntIt;

      int nUncheckedFriends = 0;
      // add friends of all surviving hits
      for(std::map<trex::TTPCAStarPoint*, float>::iterator frEl = pnt->aStarFriends.begin(); frEl != pnt->aStarFriends.end(); ++frEl){
        trex::TTPCAStarPoint* curPnt = frEl->first;

        float incCost = GetModifiedCost(frEl->second, pnt, curPnt, extendMode);
        float curCost = pnt->aStarCost + incCost;

        // veto friends that are already closed
        if (curPnt->aStarClosed) continue;
        // veto hits above current cost
        if(curCost > i){
          nUncheckedFriends ++;
          continue;
        };

        // if closer to this, add to open set and set up as child of this
        if(curPnt->aStarCost > curCost || curPnt->aStarCost < 0.){
          curPnt->aStarCost = curCost;
          curPnt->aStarParent = pnt->aStarParent;

          // only add to list if not already open
          if(!curPnt->aStarOpen){
            curPnt->aStarOpen = true;
            openSet.push_back(curPnt);
          };
        };
      };
      // if unchecked friends are still left, add this back to the open set
      if(nUncheckedFriends > 0){
        openSet.push_back(pnt);
      }
      else{
        pnt->aStarClosed = true;
      };
    };
  };

  // now add all hits associated with the given group to that group
  nGroup = 0;
  for(std::vector< trex::TTPCVolGroup >::iterator groupIt = inGroups.begin(); groupIt != inGroups.end(); ++groupIt){
    trex::TTPCVolGroup& group = *groupIt;
    trex::TTPCAStarPoint* definingPoint = definingSet[nGroup];
    ++nGroup;

    if(definingPoint != 0){
      for(std::vector<trex::TTPCAStarPoint*>::iterator pntIt = fAStarPoints.begin(); pntIt != fAStarPoints.end(); ++pntIt){
        trex::TTPCAStarPoint* pnt = *pntIt;

        // check if associated with this path
        if(pnt->aStarParent == definingPoint){
          // add vol to group
          group.AddHit(pnt->vol);
        };
      };
    };
  };
}

void trex::TTPCAStar::GetNearHitConnections(trex::TTPCVolGroupMan* volGroupMan, trex::TTPCAStarPoint* point, bool extendMode){
  // only search for friends if none are already found
  if(!point->aStarHasFriends){
    // volume group to check for nearby hits
    trex::TTPCVolGroup volCheck(fLayout);
    volCheck.AddHitMap(fHitMap);

    // find near hits and (friends)
    trex::TTPCVolGroup nearHits(fLayout);
    volGroupMan->GetNearHits(volCheck, nearHits, point->vol, trex::TTPCConnection::path);

    // save cost of connecting to each of those friends
    for(std::map<long, trex::TTPCUnitVolume*>::iterator it = nearHits.begin(); it != nearHits.end(); ++it){
      trex::TTPCAStarPoint* newPoint = 0;
      for(std::vector<trex::TTPCAStarPoint*>::iterator pnt = fAStarPoints.begin(); pnt != fAStarPoints.end(); ++pnt){
        if((*pnt)->vol == it->second){
          newPoint = *pnt;
          break;
        };
      };
      float cost = GetConnectionCost(point, newPoint);
      point->aStarFriends[newPoint] = cost;
    };
    // the current point now has friends
    point->aStarHasFriends = true;
  };
}

float trex::TTPCAStar::GetConnectionCost(trex::TTPCAStarPoint* point1, trex::TTPCAStarPoint* point2){
  // work out distances in x, y and z
  int diffX = std::abs(point1->x - point2->x);
  int diffY = std::abs(point1->y - point2->y);
  int diffZ = std::abs(point1->z - point2->z);

  // set minimum to zero
  diffX = std::max(diffX, 0);
  diffY = std::max(diffY, 0);
  diffZ = std::max(diffZ, 0);

  // work out cost from how far apart the points are in x, y and z
  float dX = (float)diffX / fXScale;
  float dY = (float)diffY / fYScale;
  float dZ = (float)diffZ / fZScale;

  float dTotal2 = dX*dX + dY*dY + dZ*dZ;
  float result = dTotal2*dTotal2;

  return result;
}
float trex::TTPCAStar::GetHeuristicCost(trex::TTPCAStarPoint* point1, trex::TTPCAStarPoint* target){
  // work out distance between point and targer in x, y and z
  float dX = (float)std::abs(point1->x - target->x) / fXScale;
  float dY = (float)std::abs(point1->y - target->y) / fYScale;
  float dZ = (float)std::abs(point1->z - target->z) / fZScale;

  return std::sqrt(dX*dX + dY*dY + dZ*dZ);
}

int trex::TTPCAStar::DoConnection(trex::TTPCAStarPoint* pathVolStart, trex::TTPCAStarPoint* pathVolEnd, bool extendMode, float maxCost){
  // return if not found
  if(!pathVolStart || !pathVolEnd) return 1;

  // define set of points to be considered at this stage in the algorithm
  std::vector<trex::TTPCAStarPoint*> openSet;
  // add the start point to the openSet

  // start point has a cumulative connection cost of 0 and is the only element in the open set
  pathVolStart->aStarCost = 0.;
  openSet.push_back(pathVolStart);

  // keep track of whether the end point has been found or not
  bool foundEnd = false;
  int nPoints = fAStarPoints.size();

  for(int i=0; i < nPoints; i++){
    // look for the next cheapest cell by minimising cost while varying id
    trex::TTPCAStarPoint* curPoint = 0;
    float curCost = 999999999.;

    if(openSet.size() == 0) break;
    for(std::vector<trex::TTPCAStarPoint*>::iterator pointIt = openSet.begin(); pointIt != openSet.end(); ++pointIt){
      trex::TTPCAStarPoint* point = *pointIt;
      // if cell has already been closed, ignore it
      if(point->aStarClosed) continue;

      // if heuristic is not defined, define it
      if(point->aStarHeuristic < 0.) point->aStarHeuristic = GetHeuristicCost(point, pathVolEnd);
      // total cost is connection cost plus heuristic cost
      float cost = point->aStarCost + point->aStarHeuristic;
      // if this is the cheapest cell, hold on to it
      if(cost < curCost){
        curPoint = point;
        curCost = cost;
      };
    };
    // if no cell was found, break out of the loop
    if(!curPoint) break;

    // if the found cell was the target, break (successfully) out of the loop
    if(curPoint == pathVolEnd){
      foundEnd = true;
      break;
    };
    // this point is now closed and erased from the open set
    curPoint->aStarClosed = true;
    openSet.erase( std::find(openSet.begin(), openSet.end(), curPoint) );

    // consider all of the current cell's friends
    for(std::map<trex::TTPCAStarPoint*, float>::iterator fr = curPoint->aStarFriends.begin(); fr != curPoint->aStarFriends.end(); ++fr){
      trex::TTPCAStarPoint* frCell = fr->first;

      if(frCell == curPoint) continue;
      if(frCell->aStarClosed) continue;

      // cost is the cost of the current cell plus the connection cost to the friend being considered
      float incCost = GetModifiedCost(fr->second, frCell, curPoint, extendMode);
      float tentativeCost = curPoint->aStarCost + incCost;

      // if this connection is the cheapest available to the friend, set the cost accordingly and set the current cell to be the friend's parent
      if(tentativeCost < frCell->aStarCost || frCell->aStarCost < 0.){
        frCell->aStarCost = tentativeCost;
        frCell->aStarParent = curPoint;
      };
      // if the current cell is already open don't bother adding it again to the list
      if(frCell->aStarOpen) continue;
      frCell->aStarOpen = true;
      openSet.push_back(frCell);

      // just break out if max cost has been hit for whatever reason
      if(maxCost > 0.){
        if(frCell->aStarCost > maxCost){
          break;
        };
      };
    };
  };

  /*int j=0;
  for(std::vector<trex::TTPCAStarPoint*>::iterator pnt = fAStarPoints.begin(); pnt != fAStarPoints.end(); ++pnt){
    ND280Debug("End point " << j << " open? " << (*pnt)->aStarOpen << ", closed? " << (*pnt)->aStarClosed);
    ND280Debug("--> getting hit connections for " << (*pnt)->vol->GetID() << ":");
    ND280Debug(" -> edge (" << (*pnt)->vol->GetEdgeX() << "," << (*pnt)->vol->GetEdgeY() << "," << (*pnt)->vol->GetEdgeZ() << ")");
    // consider all of the current cell's friends
    for(std::map<trex::TTPCAStarPoint*, int>::iterator fr = (*pnt)->aStarFriends.begin(); fr != (*pnt)->aStarFriends.end(); ++fr){
      ND280Debug(" `-> " << fr->first->vol->GetID());
    };
    j++;
    if(j>20) break;
  };
  ND280Debug("Joining " << pathVolStart->vol->GetID() << " to " << pathVolEnd->vol->GetID());*/

  // if the end wasn't found with no limit set then something has gone badly wrong
  if(!foundEnd && maxCost <= 0.){
    std::cout<<"Path not found between two points"<<std::endl;
    return 1;
  };
  return 0;
}

//MDH
//Not used
/*
std::vector< trex::THandle<trex::TTPCVolGroup> > trex::TTPCAStar::MergeGroupsAStar(std::vector< trex::THandle<trex::TTPCVolGroup> > groups, float mergeDist){
  std::vector< trex::THandle<trex::TTPCVolGroup> > inGroups = groups;
  std::vector< trex::THandle<trex::TTPCVolGroup> > mergedGroups;
  while(inGroups.size()){
    // copy id of merge from first inGroup for consistency
    unsigned int id = (*inGroups.begin())->GetID();

    std::vector< std::vector< trex::THandle<trex::TTPCVolGroup> >::iterator > mergeIts;
    mergeIts.push_back(inGroups.begin());

    // stop only when all groups overlapping are merged
    bool furtherMerges = true;
    while(furtherMerges){
      furtherMerges = false;

      // add all groups overlapping with one in the 'to merge' list
      if(inGroups.size() > 1)
      for(std::vector< std::vector< trex::THandle<trex::TTPCVolGroup> >::iterator >::iterator mergeItIt = mergeIts.begin(); mergeItIt != mergeIts.end(); ++mergeItIt){
        std::vector< trex::THandle<trex::TTPCVolGroup> >::iterator mergeIt = *mergeItIt;
        trex::THandle<trex::TTPCVolGroup> merge = *mergeIt;

        for(std::vector< trex::THandle<trex::TTPCVolGroup> >::iterator inGroupIt = inGroups.begin()+1; inGroupIt != inGroups.end(); ++inGroupIt){
          trex::THandle<trex::TTPCVolGroup> inGroup = *inGroupIt;

          bool alreadyAdded = false;
          for(std::vector< std::vector< trex::THandle<trex::TTPCVolGroup> >::iterator >::iterator mergeItIt2 = mergeIts.begin(); mergeItIt2 != mergeIts.end(); ++mergeItIt2)
          if(**mergeItIt2 == inGroup){
            alreadyAdded = true;
            break;
          };
          if(alreadyAdded) continue;

          if(FindConnectionCost(inGroup, merge, true, true, true) < mergeDist){
            mergeIts.push_back(inGroupIt);
            furtherMerges = true;
            if(furtherMerges) break;
          };
        };
        if(furtherMerges) break;
      };
    };

    // set new group with the ID of the first of the merges
    trex::THandle<trex::TTPCVolGroup> groupOut (new trex::TTPCVolGroup(fLayout, id));

    for(std::vector< std::vector< trex::THandle<trex::TTPCVolGroup> >::iterator >::iterator mergeItIt = mergeIts.begin(); mergeItIt != mergeIts.end(); ++mergeItIt){
      std::vector< trex::THandle<trex::TTPCVolGroup> >::iterator mergeIt = *mergeItIt;
      trex::THandle<trex::TTPCVolGroup> merge = *mergeIt;

      groupOut->AddHits(merge);
      groupOut->SetXLean(merge->GetXLean());
      groupOut->SetYLean(merge->GetYLean());
      groupOut->SetZLean(merge->GetZLean());

      // set to zero
      *mergeIt = trex::THandle<trex::TTPCVolGroup>();
    };

    std::vector< trex::THandle<trex::TTPCVolGroup> > newInGroups;
    for(std::vector< trex::THandle<trex::TTPCVolGroup> >::iterator grpIt = inGroups.begin(); grpIt != inGroups.end(); ++grpIt){
      trex::THandle<trex::TTPCVolGroup> grp = *grpIt;
      if(grp) newInGroups.push_back(grp);
    };
    inGroups = newInGroups;

    mergedGroups.push_back(groupOut);
  };

  return mergedGroups;
  }*/

float trex::TTPCAStar::GetModifiedCost(float cost, trex::TTPCAStarPoint* point1, trex::TTPCAStarPoint* point2, bool extendMode){
  float modifiedCost = cost;

  if(extendMode){
    modifiedCost = std::pow(cost, (float)0.25);
  };

  return modifiedCost;
}

void trex::TTPCAStar::RebootHits(){
  // reset each individual point
  for(std::vector<trex::TTPCAStarPoint*>::iterator pnt = fAStarPoints.begin(); pnt != fAStarPoints.end(); ++pnt){
    (*pnt)->aStarClosed = false;
    (*pnt)->aStarOpen = false;
    (*pnt)->aStarCost = -1.;
    (*pnt)->aStarHeuristic = -1.;
    (*pnt)->aStarParent = 0;
  };
}
