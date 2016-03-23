// eddy
#include "TTPCAStar.hxx"

ND::TTPCAStar::TTPCAStar(ND::TTPCLayout* layout){
  fLayout = layout;

  fAStarPoints = std::vector<ND::TTPCAStarPoint*>();
  fHitMap = std::map<long, ND::TTPCUnitVolume*>();

  // defualt scales and heuristic factor
  fXScale = fLayout->GetAStarXScale();
  fYScale = fLayout->GetAStarYScale();
  fZScale = fLayout->GetAStarZScale();

  fHeuristicFactor = fLayout->GetAStarHeuristicFactor();
}
ND::TTPCAStar::~TTPCAStar(){
  // delete all A* points
  for(std::vector<ND::TTPCAStarPoint*>::iterator pnt = fAStarPoints.begin(); pnt != fAStarPoints.end(); ++pnt){
    delete *pnt;
  };
}

void ND::TTPCAStar::AddHits(ND::TTPCVolGroupMan* volGroupMan, std::map<long, ND::TTPCUnitVolume*> hitMap, bool extendMode){
  // add hit map
  fHitMap = hitMap;

  // set up vector of points
  fAStarPoints = std::vector<ND::TTPCAStarPoint*>();

  // add all hits to list
  for(std::map<long, ND::TTPCUnitVolume*>::iterator hitEl = hitMap.begin(); hitEl != hitMap.end(); hitEl++){
    ND::TTPCAStarPoint* pnt = new ND::TTPCAStarPoint();

    pnt->vol = hitEl->second;

    pnt->x = hitEl->second->GetX();
    pnt->y = hitEl->second->GetY();
    pnt->z = hitEl->second->GetZ();

    pnt->aStarFriends = std::map<ND::TTPCAStarPoint*, float>();
    pnt->aStarHasFriends = false;

    fAStarPoints.push_back(pnt);
  };

  for(std::vector<ND::TTPCAStarPoint*>::iterator pointIt = fAStarPoints.begin(); pointIt != fAStarPoints.end(); ++pointIt){
    ND::TTPCAStarPoint* point = *pointIt;
    // build connections now and share friends backwards to avoid asymmetries
    GetNearHitConnections(volGroupMan, point, extendMode);
    for(std::map<ND::TTPCAStarPoint*, float>::iterator fr = point->aStarFriends.begin(); fr != point->aStarFriends.end(); ++fr){
      fr->first->aStarFriends[point] = fr->second;
    };
  };
}

void ND::TTPCAStar::AddHits(ND::TTPCAStar* prevAStar, std::map<long, ND::TTPCUnitVolume*> hitMap){
  // add hit map
  fHitMap = hitMap;

  // set up vector of points
  fAStarPoints = std::vector<ND::TTPCAStarPoint*>();

  // map to help copying points over
  std::map<ND::TTPCAStarPoint*, ND::TTPCAStarPoint*> oldNewMap;
  // copy hits that are in both the previous A* container and the map
  for(std::vector<ND::TTPCAStarPoint*>::iterator oldPntIt = prevAStar->begin(); oldPntIt != prevAStar->end(); ++oldPntIt){
    // make sure point is inside map
    ND::TTPCAStarPoint* oldPnt = *oldPntIt;
    long id = oldPnt->vol->GetID();

    if(hitMap.find(id) != hitMap.end()){
      ND::TTPCAStarPoint* pnt = new ND::TTPCAStarPoint();

      pnt->vol = oldPnt->vol;

      pnt->x = oldPnt->x;
      pnt->y = oldPnt->y;
      pnt->z = oldPnt->z;

      pnt->aStarFriends = std::map<ND::TTPCAStarPoint*, float>();

      oldNewMap[oldPnt] = pnt;
    }
  };

  // now add friends
  for(std::vector<ND::TTPCAStarPoint*>::iterator oldPntIt = prevAStar->begin(); oldPntIt != prevAStar->end(); ++oldPntIt){
    ND::TTPCAStarPoint* oldPnt = *oldPntIt;
    long id = oldPnt->vol->GetID();

    if(hitMap.find(id) != hitMap.end()){
      ND::TTPCAStarPoint* pnt = oldNewMap.at(oldPnt);
      // copy valid friends over too
      for(std::map<ND::TTPCAStarPoint*, float>::iterator fr = oldPnt->aStarFriends.begin(); fr != oldPnt->aStarFriends.end(); ++fr){
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
std::vector< ND::THandle<ND::TTPCVolGroup> > ND::TTPCAStar::ConnectVertexGroups(ND::TTPCVolGroupMan* volGroupMan, std::vector< ND::THandle<ND::TTPCVolGroup> > vertices, std::vector< ND::THandle<ND::TTPCVolGroup> > edges, int maxNo){
  std::vector< ND::THandle<ND::TTPCVolGroup> > connections;
  // find ordered sets of connections between the groups
  std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > orderedConnections = ConnectVertexGroupsOrdered(volGroupMan, vertices, edges, maxNo);

  // convert the ordered groups to unordered groups
  for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator path = orderedConnections.begin(); path != orderedConnections.end(); ++path){
    ND::THandle<ND::TTPCVolGroup> connection(new ND::TTPCVolGroup(fLayout));
    for(std::vector<ND::TTPCPathVolume*>::iterator cell = (*path)->begin(); cell != (*path)->end(); ++cell){
      connection->AddHit((*cell)->GetUnitVolume());
    };
    connections.push_back(connection);
  };

  return connections;
}
*/

void ND::TTPCAStar::ConnectVertexGroupsOrdered(ND::TTPCVolGroupMan* volGroupMan, std::vector< ND::TTPCVolGroup >& vertices, std::vector< ND::TTPCVolGroup >& edges, std::vector< ND::TTPCOrderedVolGroup >& connections, int maxNo){

  // keep track of number of iterations to break when exceeded
  int i=0;
  // loop over each vertex and each edge
  for(std::vector< ND::TTPCVolGroup >::iterator ver = vertices.begin(); ver != vertices.end(); ++ver){
    for(std::vector< ND::TTPCVolGroup >::iterator grp = edges.begin(); grp != edges.end(); ++grp){
      // add connection between the groups to set of connections
      connections.push_back(fLayout);
      ConnectGroupPair(*ver,*grp,connections.back(), true,false);
      if(connections.back().empty()) connections.pop_back();
      // break out if the number of iterations is too high
      i++;
      if(i>=maxNo) break;
    };
    if(i>=maxNo) break;
  };

  // connect vertices to each other as well
  std::vector< ND::TTPCOrderedVolGroup > vertexConnections;
  ConnectGroupsOrdered(volGroupMan, vertices, vertexConnections, true, false, maxNo);
  // add inter-vertex connections to group
  for(std::vector< ND::TTPCOrderedVolGroup >::iterator it = vertexConnections.begin(); it != vertexConnections.end(); ++it){
    connections.push_back(*it);
  };

}

//MDH
//Not used
/*
std::vector< ND::THandle<ND::TTPCVolGroup> > ND::TTPCAStar::ConnectGroups(ND::TTPCVolGroupMan* volGroupMan, std::vector< ND::THandle<ND::TTPCVolGroup> > groups, bool allConnections, int maxNo){
  std::vector< ND::THandle<ND::TTPCVolGroup> > connections;
  // find ordered sets of connections between the groups
  std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > orderedConnections = ConnectGroupsOrdered(volGroupMan, groups, false, allConnections, maxNo);

  // convert the ordered groups to unordered groups
  for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator path = orderedConnections.begin(); path != orderedConnections.end(); ++path){
    ND::THandle<ND::TTPCVolGroup> connection(new ND::TTPCVolGroup(fLayout));
    for(std::vector<ND::TTPCPathVolume*>::iterator cell = (*path)->begin(); cell != (*path)->end(); ++cell){
      connection->AddHit((*cell)->GetUnitVolume());
    };
    connections.push_back(connection);
  };

  return connections;
}
*/

void ND::TTPCAStar::ConnectGroupsOrdered(ND::TTPCVolGroupMan* volGroupMan, std::vector< ND::TTPCVolGroup >& groups, std::vector< ND::TTPCOrderedVolGroup >& connections, bool vertices, bool allConnections, int maxNo){

  // keep track of number of iterations to break when exceeded
  int i=0;
  // loop over each pair of groups once
  for(std::vector< ND::TTPCVolGroup >::iterator grp1 = groups.begin(); grp1 != groups.end(); ++grp1){
    std::vector< ND::TTPCVolGroup >::iterator grp2Begin = allConnections ? groups.begin() : grp1+1;
    std::vector< ND::TTPCVolGroup >::iterator grp2End = groups.end();
    for(std::vector< ND::TTPCVolGroup >::iterator grp2 = grp2Begin; grp2 != grp2End; ++grp2){
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
void ND::TTPCAStar::ClearRedundancies(ND::TTPCVolGroupMan* volGroupMan, std::vector< ND::TTPCVolGroup >& groups, int maxNo){
  // the next block of code is pretty ugly but needed to remove bugs in some very specific topologies
  // first, merge groups which are too close together

  std::vector< ND::TTPCVolGroup > mergedGroups;

  while(groups.size()){
    // copy id of merge from first inGroup for consistency
    unsigned int id = (*groups.begin()).GetID();

    std::vector< std::vector< ND::TTPCVolGroup >::iterator > mergeIts;
    mergeIts.push_back(groups.begin());

    // stop only when all groups overlapping are merged
    bool furtherMerges = true;
    while(furtherMerges){
      furtherMerges = false;

      // add all groups overlapping with one in the 'to merge' list
      if(groups.size() > 1)
      for(std::vector< std::vector< ND::TTPCVolGroup >::iterator >::iterator mergeItIt = mergeIts.begin(); mergeItIt != mergeIts.end(); ++mergeItIt){
        std::vector< ND::TTPCVolGroup >::iterator mergeIt = *mergeItIt;
        ND::TTPCVolGroup& merge = *mergeIt;

        for(std::vector< ND::TTPCVolGroup >::iterator inGroupIt = groups.begin(); inGroupIt != groups.end(); ++inGroupIt){
          ND::TTPCVolGroup& inGroup = *inGroupIt;

          bool alreadyAdded = false;
          for(std::vector< std::vector< ND::TTPCVolGroup >::iterator >::iterator mergeItIt2 = mergeIts.begin(); mergeItIt2 != mergeIts.end(); ++mergeItIt2)
          if(*mergeItIt2 == inGroupIt){
            alreadyAdded = true;
            break;
          };
          if(alreadyAdded) continue;

          if( volGroupMan->GetGroupGroupOverlap(inGroup, merge, ND::TTPCConnection::edgeMerge, true) ){
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
    ND::TTPCVolGroup& groupOut = mergedGroups.back();
    std::vector< std::vector<ND::TTPCVolGroup>::iterator > delGroups;

    for(std::vector< std::vector< ND::TTPCVolGroup >::iterator >::iterator mergeItIt = mergeIts.begin(); mergeItIt != mergeIts.end(); ++mergeItIt){
      std::vector< ND::TTPCVolGroup >::iterator mergeIt = *mergeItIt;
      ND::TTPCVolGroup& merge = *mergeIt;

      groupOut.AddHits(merge);
      groupOut.SetXLean(merge.GetXLean());
      groupOut.SetYLean(merge.GetYLean());
      groupOut.SetZLean(merge.GetZLean());

      delGroups.push_back(mergeIt);
    };

    std::vector< ND::TTPCVolGroup > newInGroups;
    for(std::vector< ND::TTPCVolGroup >::iterator grpIt = groups.begin(); grpIt != groups.end(); ++grpIt){
      if(std::find(delGroups.begin(),delGroups.end(),grpIt)==delGroups.end()){
      newInGroups.push_back(*grpIt);
    }
      groups = std::move(newInGroups);
    }
  }

  // parallel arrays for groups and list of groups they overlap
  std::vector< std::vector<int> > connectedGroups;
  for(std::vector< ND::TTPCVolGroup >::iterator grpIt = mergedGroups.begin(); grpIt != mergedGroups.end(); ++grpIt){
    connectedGroups.emplace_back();
  }
  std::vector<bool> groupsMaySurvive (mergedGroups.size(), true);

  // keep track of number of iterations to break when exceeded
  int nGroups = (int)mergedGroups.size();
  int count=0;

  // loop over each pair of groups once
  for(int i=0; i<nGroups; ++i){
    for(int j=i+1; j<nGroups; ++j){
      ND::TTPCVolGroup& grp1 = mergedGroups[i];
      ND::TTPCVolGroup& grp2 = mergedGroups[j];
      // add connection between the groups to set of connections
      ND::TTPCOrderedVolGroup path(fLayout);
      ConnectGroupPair(grp1, grp2,path);
      // kill any groups associated with this path other than its direct ends
      for(int k=0; k<nGroups; ++k){
        if(k == i) continue;
        if(k == j) continue;
        ND::TTPCVolGroup& grp3 = mergedGroups[k];
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
  std::vector< ND::TTPCVolGroup > survivors;

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

  groups=std::move(survivors);
}

void ND::TTPCAStar::ClearVertexConnectionRedundancies(ND::TTPCVolGroupMan* volGroupMan, std::vector< ND::TTPCOrderedVolGroup >& paths, std::vector< ND::TTPCVolGroup >& vertices){
  std::vector< ND::TTPCOrderedVolGroup > outPaths;

  for(std::vector< ND::TTPCOrderedVolGroup >::iterator path = paths.begin(); path != paths.end(); ++path){
    bool redundancy = false;
    for(std::vector< ND::TTPCVolGroup >::iterator vert = vertices.begin(); vert != vertices.end(); ++vert){
      if(volGroupMan->GetPathVolOverlap(*path, vert->GetAverageVol(), ND::TTPCConnection::vertexPath)){
	redundancy=true;
        break;
      }
    }
    if(!redundancy){
      outPaths.push_back(*path);
    }
  }

  paths=std::move(outPaths);
  
}

void ND::TTPCAStar::ConnectGroupPair(ND::TTPCVolGroup& group1, ND::TTPCVolGroup& group2, ND::TTPCOrderedVolGroup& connection, bool vertexGroup1, bool vertexGroup2, bool fullASICPenalty, bool extendMode){

  // get unique ids and A* indices
  ND::TTPCUnitVolume* startCell = group1.GetAverageVol();
  ND::TTPCUnitVolume* endCell = group2.GetAverageVol();

  if(!startCell || !endCell) return;

  // reset relevant variables each time a connection needs to be made
  RebootHits();

  // assign A* points based on these
  ND::TTPCAStarPoint* startPoint=0;
  ND::TTPCAStarPoint* endPoint=0;
  for(std::vector<ND::TTPCAStarPoint*>::iterator pnt = fAStarPoints.begin(); pnt != fAStarPoints.end(); ++pnt){
    if((*pnt)->vol == startCell) startPoint = *pnt;
    if((*pnt)->vol == endCell) endPoint = *pnt;
  };

  // set up the A* variables chaining start point to end point
  int connectionReturnCode = DoConnection(startPoint, endPoint, fullASICPenalty, extendMode);
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
  ND::TTPCAStarPoint* curPoint = endPoint;
  int nPoints = fAStarPoints.size();
  for(int i=0; i < nPoints; i++){
    connection.AddCell(curPoint->vol);

    ND::TTPCAStarPoint* nextPoint = curPoint->aStarParent;
    if(!nextPoint) break; 
    curPoint = nextPoint;
  };
  // add the first cell so long as the connection isn't one cell long
  if(curPoint != endPoint) connection.AddCell(curPoint->vol);

}

float ND::TTPCAStar::FindConnectionCost(ND::TTPCVolGroup& group1, ND::TTPCVolGroup& group2, bool fullASICPenalty, bool extendMode, bool reduced, float maxCost){
  ND::TTPCUnitVolume* vol1 = group1.GetAverageVol();
  ND::TTPCUnitVolume* vol2 = group2.GetAverageVol();

  return FindConnectionCost(vol1, vol2, fullASICPenalty, extendMode, reduced, maxCost);
}
float ND::TTPCAStar::FindConnectionCost(ND::TTPCVolGroup& group, ND::TTPCUnitVolume* vol, bool fullASICPenalty, bool extendMode, bool reduced, float maxCost){
  ND::TTPCUnitVolume* groupVol = group.GetAverageVol();

  return FindConnectionCost(groupVol, vol, fullASICPenalty, extendMode, reduced, maxCost);
}
float ND::TTPCAStar::FindConnectionCost(ND::TTPCUnitVolume* vol1, ND::TTPCUnitVolume* vol2, bool fullASICPenalty, bool extendMode, bool reduced, float maxCost){
  float cost = -1.;
  if(!vol1 || !vol2) return cost;

  if(reduced){
    maxCost *= fLayout->GetAStarPathologyPenalty();
  }

  // reset relevant variables each time a connection needs to be made
  RebootHits();

  // assign A* points based on these
  ND::TTPCAStarPoint* startPoint=0;
  ND::TTPCAStarPoint* endPoint=0;
  for(std::vector<ND::TTPCAStarPoint*>::iterator pnt = fAStarPoints.begin(); pnt != fAStarPoints.end(); ++pnt){
    if((*pnt)->vol == vol1) startPoint = *pnt;
    if((*pnt)->vol == vol2) endPoint = *pnt;
  };

  int connectionReturnCode = DoConnection(startPoint, endPoint, fullASICPenalty, extendMode, maxCost);
  if(connectionReturnCode){
    return cost;
  }
  else{
    cost = endPoint->aStarCost;

    if(reduced){
      // reduce cost depending on number of pathological hits in the path
      int nHits = 0;
      int nPHits = 0;
      int nNHits = 0;
      // start at end index and add the chain of cells connecting it from the start to the connection
      ND::TTPCAStarPoint* curPoint = endPoint;
      int nPoints = fAStarPoints.size();
      for(int i=0; i < nPoints; i++){
        nHits++;
        if(curPoint->vol->GetPathology()){
          nPHits++;
        }
        else{
          nNHits++;
        };

        ND::TTPCAStarPoint* nextPoint = curPoint->aStarParent;
        if(!nextPoint) break;
        curPoint = nextPoint;
      };

      if(nHits > 0){
        float reduceFactor = ((float)nNHits/(float)nHits) + (((float)nPHits/(float)nHits) / fLayout->GetAStarPathologyPenalty());
        cost *= reduceFactor;
      };
    };
    if(cost <= 0.){
      cost = maxCost;
    };
  };

  return cost;
}

void ND::TTPCAStar::AssociateBestHits(ND::TTPCVolGroupMan* volGroupMan, std::vector< ND::TTPCOrderedVolGroup >& inPaths, bool fullASICPenalty, bool extendMode, float maxCost){
  // add hits in path to group, making sure the same isn't added twice
  std::vector< ND::TTPCVolGroup > hitsGroups;
  std::set<ND::TTPCUnitVolume*> volsAdded;
  for(std::vector< ND::TTPCOrderedVolGroup >::iterator inPathIt = inPaths.begin(); inPathIt != inPaths.end(); ++inPathIt){
    ND::TTPCOrderedVolGroup& inPath = *inPathIt;

    hitsGroups.emplace_back(fLayout);
    ND::TTPCVolGroup& hitsGroup = hitsGroups.back();
    for(std::vector<ND::TTPCPathVolume*>::iterator pathVolIt = inPath.begin(); pathVolIt != inPath.end(); ++pathVolIt){
      ND::TTPCPathVolume* pathVol = *pathVolIt;
      ND::TTPCUnitVolume* vol = pathVol->GetUnitVolume();

      if(!volsAdded.count(vol)){
        hitsGroup.AddHit(vol);
        volsAdded.insert(vol);
      };
    };

  };

  // merge best hits in groups
  MergeBestHits(volGroupMan, hitsGroups, fullASICPenalty, extendMode, maxCost);

  // add merged hits to appropriate paths
  for(unsigned int i=0; i<inPaths.size(); ++i){
    ND::TTPCOrderedVolGroup& inPath = inPaths.at(i);
    ND::TTPCVolGroup& hitsGroup = hitsGroups.at(i);

    inPath.AddExtendedHits(hitsGroup);
  };
}

//MDH
//Not used
/*
void ND::TTPCAStar::MergeBestHits(ND::TTPCVolGroupMan* volGroupMan, ND::THandle<ND::TTPCVolGroup> inGroup, bool fullASICPenalty, bool extendMode, float maxCost){
  std::vector< ND::THandle<ND::TTPCVolGroup> > dummyVector;
  dummyVector.push_back(inGroup);

  MergeBestHits(volGroupMan, dummyVector, fullASICPenalty, extendMode, maxCost);
}
*/

void ND::TTPCAStar::MergeBestHits(ND::TTPCVolGroupMan* volGroupMan, std::vector< ND::TTPCVolGroup >& inGroups, bool fullASICPenalty, bool extendMode, float maxCost){
  // reset hits
  RebootHits();

  // define container of all near hits currently under consideration, starting with path hits
  std::vector<ND::TTPCAStarPoint*> openSet;
  std::vector<ND::TTPCAStarPoint*> definingSet (inGroups.size(), 0);
  int nGroup = 0;
  for(std::vector< ND::TTPCVolGroup >::iterator groupIt = inGroups.begin(); groupIt != inGroups.end(); ++groupIt){
    ND::TTPCVolGroup& group = *groupIt;

    for(std::map<long, ND::TTPCUnitVolume*>::iterator hitEl = group.begin(); hitEl != group.end(); ++hitEl){
      ND::TTPCUnitVolume* vol = hitEl->second;

      for(std::vector<ND::TTPCAStarPoint*>::iterator pntIt = fAStarPoints.begin(); pntIt != fAStarPoints.end(); ++pntIt){
        ND::TTPCAStarPoint* pnt = *pntIt;

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
    std::vector<ND::TTPCAStarPoint*> curOpenSet;
    for(std::vector<ND::TTPCAStarPoint*>::iterator pntIt = openSet.begin(); pntIt != openSet.end(); ++pntIt){
      ND::TTPCAStarPoint* pnt = *pntIt;

      if(!pnt->aStarClosed){
        curOpenSet.push_back(*pntIt);
      };
    };
    openSet.clear();

    // loop once over current set
    for(std::vector<ND::TTPCAStarPoint*>::iterator pntIt = curOpenSet.begin(); pntIt != curOpenSet.end(); ++pntIt){
      ND::TTPCAStarPoint* pnt = *pntIt;

      int nUncheckedFriends = 0;
      // add friends of all surviving hits
      for(std::map<ND::TTPCAStarPoint*, float>::iterator frEl = pnt->aStarFriends.begin(); frEl != pnt->aStarFriends.end(); ++frEl){
        ND::TTPCAStarPoint* curPnt = frEl->first;

        float incCost = GetModifiedCost(frEl->second, pnt, curPnt, fullASICPenalty, extendMode);
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
  for(std::vector< ND::TTPCVolGroup >::iterator groupIt = inGroups.begin(); groupIt != inGroups.end(); ++groupIt){
    ND::TTPCVolGroup& group = *groupIt;
    ND::TTPCAStarPoint* definingPoint = definingSet[nGroup];
    ++nGroup;

    if(definingPoint != 0){
      for(std::vector<ND::TTPCAStarPoint*>::iterator pntIt = fAStarPoints.begin(); pntIt != fAStarPoints.end(); ++pntIt){
        ND::TTPCAStarPoint* pnt = *pntIt;

        // check if associated with this path
        if(pnt->aStarParent == definingPoint){
          // add vol to group
          group.AddHit(pnt->vol);
        };
      };
    };
  };
}

void ND::TTPCAStar::GetNearHitConnections(ND::TTPCVolGroupMan* volGroupMan, ND::TTPCAStarPoint* point, bool extendMode){
  // only search for friends if none are already found
  if(!point->aStarHasFriends){
    // volume group to check for nearby hits
    ND::TTPCVolGroup volCheck(fLayout);
    volCheck.AddHitMap(fHitMap);

    // find near hits and (friends)
    ND::TTPCVolGroup nearHits(fLayout);
    volGroupMan->GetNearHits(volCheck, nearHits, point->vol, ND::TTPCConnection::path);

    // save cost of connecting to each of those friends
    for(std::map<long, ND::TTPCUnitVolume*>::iterator it = nearHits.begin(); it != nearHits.end(); ++it){
      ND::TTPCAStarPoint* newPoint = 0;
      for(std::vector<ND::TTPCAStarPoint*>::iterator pnt = fAStarPoints.begin(); pnt != fAStarPoints.end(); ++pnt){
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

float ND::TTPCAStar::GetConnectionCost(ND::TTPCAStarPoint* point1, ND::TTPCAStarPoint* point2, bool penalty){
  // work out distances in x, y and z
  int diffX = std::abs(point1->x - point2->x);
  int diffY = std::abs(point1->y - point2->y);
  int diffZ = std::abs(point1->z - point2->z);

  // reduce costs if jumping gaps
  if(point1->vol->GetSegX() != point2->vol->GetSegX()){
    diffX -= fLayout->GetGapOffsetX();
  };
  if(point1->vol->GetSegY() != point2->vol->GetSegY()){
    diffY -= fLayout->GetGapOffsetY();
  };
  if(point1->vol->GetSegZ() != point2->vol->GetSegZ()){
    diffZ -= fLayout->GetGapOffsetZ();
  };

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

  // if demmanding a low charge penalty or a pathology penalty, add it in
  if(penalty){
    int nPathology = (int)point1->vol->GetPathology() + (int)point2->vol->GetPathology();
    float penaltyVal = fLayout->GetAStarPathologyPenalty() * (float)nPathology;
    result *= penaltyVal;
  };

  return result;
}
float ND::TTPCAStar::GetHeuristicCost(ND::TTPCAStarPoint* point1, ND::TTPCAStarPoint* target){
  // work out distance between point and targer in x, y and z
  float dX = (float)std::abs(point1->x - target->x) / fXScale;
  float dY = (float)std::abs(point1->y - target->y) / fYScale;
  float dZ = (float)std::abs(point1->z - target->z) / fZScale;

  return std::sqrt(dX*dX + dY*dY + dZ*dZ);
}

int ND::TTPCAStar::DoConnection(ND::TTPCAStarPoint* pathVolStart, ND::TTPCAStarPoint* pathVolEnd, bool fullASICPenalty, bool extendMode, float maxCost){
  // return if not found
  if(!pathVolStart || !pathVolEnd) return 1;

  // define set of points to be considered at this stage in the algorithm
  std::vector<ND::TTPCAStarPoint*> openSet;
  // add the start point to the openSet

  // start point has a cumulative connection cost of 0 and is the only element in the open set
  pathVolStart->aStarCost = 0.;
  openSet.push_back(pathVolStart);

  // keep track of whether the end point has been found or not
  bool foundEnd = false;
  int nPoints = fAStarPoints.size();

  for(int i=0; i < nPoints; i++){
    // look for the next cheapest cell by minimising cost while varying id
    ND::TTPCAStarPoint* curPoint = 0;
    float curCost = 999999999.;

    if(openSet.size() == 0) break;
    for(std::vector<ND::TTPCAStarPoint*>::iterator pointIt = openSet.begin(); pointIt != openSet.end(); ++pointIt){
      ND::TTPCAStarPoint* point = *pointIt;
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
    for(std::map<ND::TTPCAStarPoint*, float>::iterator fr = curPoint->aStarFriends.begin(); fr != curPoint->aStarFriends.end(); ++fr){
      ND::TTPCAStarPoint* frCell = fr->first;

      if(frCell == curPoint) continue;
      if(frCell->aStarClosed) continue;

      // cost is the cost of the current cell plus the connection cost to the friend being considered
      float incCost = GetModifiedCost(fr->second, frCell, curPoint, fullASICPenalty, extendMode);
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
  for(std::vector<ND::TTPCAStarPoint*>::iterator pnt = fAStarPoints.begin(); pnt != fAStarPoints.end(); ++pnt){
    ND280Debug("End point " << j << " open? " << (*pnt)->aStarOpen << ", closed? " << (*pnt)->aStarClosed);
    ND280Debug("--> getting hit connections for " << (*pnt)->vol->GetID() << ":");
    ND280Debug(" -> edge (" << (*pnt)->vol->GetEdgeX() << "," << (*pnt)->vol->GetEdgeY() << "," << (*pnt)->vol->GetEdgeZ() << ")");
    // consider all of the current cell's friends
    for(std::map<ND::TTPCAStarPoint*, int>::iterator fr = (*pnt)->aStarFriends.begin(); fr != (*pnt)->aStarFriends.end(); ++fr){
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
std::vector< ND::THandle<ND::TTPCVolGroup> > ND::TTPCAStar::MergeGroupsAStar(std::vector< ND::THandle<ND::TTPCVolGroup> > groups, float mergeDist){
  std::vector< ND::THandle<ND::TTPCVolGroup> > inGroups = groups;
  std::vector< ND::THandle<ND::TTPCVolGroup> > mergedGroups;
  while(inGroups.size()){
    // copy id of merge from first inGroup for consistency
    unsigned int id = (*inGroups.begin())->GetID();

    std::vector< std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator > mergeIts;
    mergeIts.push_back(inGroups.begin());

    // stop only when all groups overlapping are merged
    bool furtherMerges = true;
    while(furtherMerges){
      furtherMerges = false;

      // add all groups overlapping with one in the 'to merge' list
      if(inGroups.size() > 1)
      for(std::vector< std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator >::iterator mergeItIt = mergeIts.begin(); mergeItIt != mergeIts.end(); ++mergeItIt){
        std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator mergeIt = *mergeItIt;
        ND::THandle<ND::TTPCVolGroup> merge = *mergeIt;

        for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator inGroupIt = inGroups.begin()+1; inGroupIt != inGroups.end(); ++inGroupIt){
          ND::THandle<ND::TTPCVolGroup> inGroup = *inGroupIt;

          bool alreadyAdded = false;
          for(std::vector< std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator >::iterator mergeItIt2 = mergeIts.begin(); mergeItIt2 != mergeIts.end(); ++mergeItIt2)
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
    ND::THandle<ND::TTPCVolGroup> groupOut (new ND::TTPCVolGroup(fLayout, id));

    for(std::vector< std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator >::iterator mergeItIt = mergeIts.begin(); mergeItIt != mergeIts.end(); ++mergeItIt){
      std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator mergeIt = *mergeItIt;
      ND::THandle<ND::TTPCVolGroup> merge = *mergeIt;

      groupOut->AddHits(merge);
      groupOut->SetXLean(merge->GetXLean());
      groupOut->SetYLean(merge->GetYLean());
      groupOut->SetZLean(merge->GetZLean());

      // set to zero
      *mergeIt = ND::THandle<ND::TTPCVolGroup>();
    };

    std::vector< ND::THandle<ND::TTPCVolGroup> > newInGroups;
    for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator grpIt = inGroups.begin(); grpIt != inGroups.end(); ++grpIt){
      ND::THandle<ND::TTPCVolGroup> grp = *grpIt;
      if(grp) newInGroups.push_back(grp);
    };
    inGroups = newInGroups;

    mergedGroups.push_back(groupOut);
  };

  return mergedGroups;
  }*/

float ND::TTPCAStar::GetModifiedCost(float cost, ND::TTPCAStarPoint* point1, ND::TTPCAStarPoint* point2, bool fullASICPenalty, bool extendMode){
  float penaltyVal = fLayout->GetAStarAssociatePathologyPenalty();
  float modifiedCost = cost;

  if(extendMode){
    modifiedCost = std::pow(cost, (float)0.25);
  };
  if(fullASICPenalty){
    if(point1->vol->GetFullASICTagged() && point2->vol->GetFullASICTagged()){
      modifiedCost *= penaltyVal;
    };
  };

  return modifiedCost;
}

void ND::TTPCAStar::RebootHits(){
  // reset each individual point
  for(std::vector<ND::TTPCAStarPoint*>::iterator pnt = fAStarPoints.begin(); pnt != fAStarPoints.end(); ++pnt){
    (*pnt)->aStarClosed = false;
    (*pnt)->aStarOpen = false;
    (*pnt)->aStarCost = -1.;
    (*pnt)->aStarHeuristic = -1.;
    (*pnt)->aStarParent = 0;
  };
}
