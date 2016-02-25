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
std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > ND::TTPCAStar::ConnectVertexGroupsOrdered(ND::TTPCVolGroupMan* volGroupMan, std::vector< ND::THandle<ND::TTPCVolGroup> > vertices, std::vector< ND::THandle<ND::TTPCVolGroup> > edges, int maxNo){
  std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > connections;

  // keep track of number of iterations to break when exceeded
  int i=0;
  // loop over each vertex and each edge
  for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator ver = vertices.begin(); ver != vertices.end(); ++ver){
    for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator grp = edges.begin(); grp != edges.end(); ++grp){
      // add connection between the groups to set of connections
      ND::THandle<ND::TTPCOrderedVolGroup> connection = ConnectGroupPair(*ver,*grp, true,false);
      if(!connection->empty()) connections.push_back(connection);
      // break out if the number of iterations is too high
      i++;
      if(i>=maxNo) break;
    };
    if(i>=maxNo) break;
  };

  // connect vertices to each other as well
  std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > vertexConnections = ConnectGroupsOrdered(volGroupMan, vertices, true, maxNo);
  // add inter-vertex connections to group
  for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator it = vertexConnections.begin(); it != vertexConnections.end(); ++it){
    connections.push_back(*it);
  };

  return connections;
}

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
std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > ND::TTPCAStar::ConnectGroupsOrdered(ND::TTPCVolGroupMan* volGroupMan, std::vector< ND::THandle<ND::TTPCVolGroup> > groups, bool vertices, bool allConnections, int maxNo){
  std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > connections;

  // keep track of number of iterations to break when exceeded
  int i=0;
  // loop over each pair of groups once
  for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator grp1 = groups.begin(); grp1 != groups.end(); ++grp1){
    std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator grp2Begin = allConnections ? groups.begin() : grp1+1;
    std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator grp2End = groups.end();
    for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator grp2 = grp2Begin; grp2 != grp2End; ++grp2){
      if(*grp1 == *grp2) continue;
      i++;
      // add connection between the groups to set of connections
      ND::THandle<ND::TTPCOrderedVolGroup> connection = ConnectGroupPair(*grp1,*grp2, vertices,vertices);
      if(!connection->empty()) connections.push_back(connection);
      // break out if the number of iterations is too high
      if(i>=maxNo) break;
    };
    if(i>=maxNo) break;
  };

  return connections;
}
std::vector< ND::THandle<ND::TTPCVolGroup> > ND::TTPCAStar::ClearRedundancies(ND::TTPCVolGroupMan* volGroupMan, std::vector< ND::THandle<ND::TTPCVolGroup> > groups, int maxNo){
  // the next block of code is pretty ugly but needed to remove bugs in some very specific topologies
  // first, merge groups which are too close together
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

        for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator inGroupIt = inGroups.begin(); inGroupIt != inGroups.end(); ++inGroupIt){
          ND::THandle<ND::TTPCVolGroup> inGroup = *inGroupIt;

          bool alreadyAdded = false;
          for(std::vector< std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator >::iterator mergeItIt2 = mergeIts.begin(); mergeItIt2 != mergeIts.end(); ++mergeItIt2)
          if(**mergeItIt2 == inGroup){
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

  // parallel arrays for groups and list of groups they overlap
  std::vector< ND::THandle<ND::TTPCVolGroup> > findGroups (mergedGroups);
  std::vector< std::vector<int> > connectedGroups;
  for(std::vector<ND::THandle<ND::TTPCVolGroup> >::iterator grpIt = findGroups.begin(); grpIt != findGroups.end(); ++grpIt) connectedGroups.push_back(std::vector<int>());
  std::vector<bool> groupsMaySurvive (findGroups.size(), true);

  // keep track of number of iterations to break when exceeded
  int nGroups = (int)findGroups.size();
  int count=0;

  // loop over each pair of groups once
  for(int i=0; i<nGroups; ++i){
    for(int j=i+1; j<nGroups; ++j){
      ND::THandle<ND::TTPCVolGroup> grp1 = findGroups[i];
      ND::THandle<ND::TTPCVolGroup> grp2 = findGroups[j];
      // add connection between the groups to set of connections
      ND::THandle<ND::TTPCOrderedVolGroup> path = ConnectGroupPair(grp1, grp2);
      // kill any groups associated with this path other than its direct ends
      for(int k=0; k<nGroups; ++k){
        if(k == i) continue;
        if(k == j) continue;
        ND::THandle<ND::TTPCVolGroup> grp3 = findGroups[k];
        if(volGroupMan->GetPathVolOverlap(path, grp3->GetAverageVol())){
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
  std::vector< ND::THandle<ND::TTPCVolGroup> > survivors = std::vector< ND::THandle<ND::TTPCVolGroup> >();

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
      survivors.push_back(findGroups[i]);

      // mark everything this group overlapped with for death
      int connection2Size = connectedGroups[i].size();
      for(int conn=0; conn<connection2Size; ++conn) groupsMaySurvive[ connectedGroups[i][conn] ] = false;
      connectedGroups[i] = std::vector<int>();
    };
  };

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
            ND::THandle<ND::TTPCVolGroup> grp = volGroupMan->MergeGroups(findGroups[merge1], findGroups[merge2]);
            survivors.push_back(grp);
          };
          break;
        };
      };
    };
  };

  return survivors;
}
std::vector< ND::THandle<ND::TTPCVolGroup> > ND::TTPCAStar::ExperimentalClearRedundancies(ND::TTPCVolGroupMan* volGroupMan, std::vector< ND::THandle<ND::TTPCVolGroup> > groups, bool forceReconstructable, int maxNo){
  // the next block of code is pretty ugly but needed to remove bugs in some very specific topologies
  // merge groups which are too close together
  int edgeClearDist = fLayout->GetAltEdgeHitConnectDist();

  // parallel arrays for groups and list of groups they overlap
  std::vector< ND::THandle<ND::TTPCVolGroup> > mergedGroups = MergeGroupsAStar(groups, edgeClearDist);
  //std::vector< ND::THandle<ND::TTPCVolGroup> > mergedGroups = groups;
  std::vector< std::vector<int> > connectedGroups;

  // if forcing a reconstructable output, require at least two ends
  if(forceReconstructable){
    if(mergedGroups.size() < 2){
      std::vector< ND::THandle<ND::TTPCVolGroup> > survivors;

      // add most separated pair
      float maxSep = -1.;
      ND::THandle<ND::TTPCVolGroup> grp1;
      ND::THandle<ND::TTPCVolGroup> grp2;

      if(groups.size() > 1){
        // loop over each pair of groups once
        for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator grp1It = groups.begin(); grp1It != groups.end(); ++grp1It){
          ND::THandle<ND::TTPCVolGroup> tempGrp1 = *grp1It;
          if(!tempGrp1) continue;

          for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator grp2It = grp1It+1; grp2It != groups.end(); ++grp2It){
            ND::THandle<ND::TTPCVolGroup> tempGrp2 = *grp2It;
            if(!tempGrp2) continue;

            float sep = FindConnectionCost(tempGrp1, tempGrp2, true, true, true);
            if(sep > maxSep){
              sep = maxSep;
              grp1 = tempGrp1;
              grp2 = tempGrp2;
            };
          };
        };
      };

      // add most separated
      if(maxSep >= 0.){
        survivors.push_back(grp1);
        survivors.push_back(grp2);
      };

      return survivors;
    };
  };

  for(std::vector<ND::THandle<ND::TTPCVolGroup> >::iterator grpIt = mergedGroups.begin(); grpIt != mergedGroups.end(); ++grpIt) connectedGroups.push_back(std::vector<int>());
  std::vector<bool> groupsMaySurvive (mergedGroups.size(), true);

  // keep track of number of iterations to break when exceeded
  int nGroups = (int)mergedGroups.size();
  int count=0;

  // loop over each pair of groups once
  for(int i=0; i<nGroups; ++i){
    for(int j=i+1; j<nGroups; ++j){
      ND::THandle<ND::TTPCVolGroup> grp1 = mergedGroups[i];
      ND::THandle<ND::TTPCVolGroup> grp2 = mergedGroups[j];
      // add connection between the groups to set of connections
      ND::THandle<ND::TTPCOrderedVolGroup> path = ConnectGroupPair(grp1, grp2);
      // kill any groups associated with this path other than its direct ends
      for(int k=0; k<nGroups; ++k){
        if(k == i) continue;
        if(k == j) continue;
        ND::THandle<ND::TTPCVolGroup> grp3 = mergedGroups[k];
        bool found = false;

        for(std::vector<ND::TTPCPathVolume*>::iterator pathVolIt = path->begin(); pathVolIt != path->end(); ++pathVolIt){
          ND::TTPCPathVolume* pathVol = *pathVolIt;
          ND::TTPCUnitVolume* vol = pathVol->GetUnitVolume();

          if(FindConnectionCost(grp3, vol, true, true, true) <= edgeClearDist){
            found = true;
            break;
          };
        };
        if(found){
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
  std::vector< ND::THandle<ND::TTPCVolGroup> > survivors = std::vector< ND::THandle<ND::TTPCVolGroup> >();

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
      connectedGroups[i] = std::vector<int>();
    };
  };

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
            ND::THandle<ND::TTPCVolGroup> grp = volGroupMan->MergeGroups(mergedGroups[merge1], mergedGroups[merge2]);
            survivors.push_back(grp);
          };
          break;
        };
      };
    };
  };

  return survivors;
}
std::vector< ND::THandle<ND::TTPCVolGroup> > ND::TTPCAStar::ExperimentalClearRedundancies2(ND::TTPCVolGroupMan* volGroupMan, std::vector< ND::THandle<ND::TTPCVolGroup> > groups, int maxNo){
  // don't check if there are only two
  if(groups.size() < 3){
    return groups;
  };

  // merge groups which are too close together
  int edgeClearDist = fLayout->GetAltEdgeHitConnectDist();

  std::vector< ND::THandle<ND::TTPCVolGroup> > outGroups;

  // make list of all cross redundancies, and summed redundancy for each group
  int nGroups = (int)groups.size();
  std::map< std::pair<int, int>, std::vector<int> > crossRedundancies;
  /*std::vector< std::vector<int> > conservativeCrossRedundancies;
  for(int i=0; i<nGroups; ++i){
    conservativeCrossRedundancies.push_back( std::vector<int>() );
  };*/

  for(int i=0; i<nGroups; ++i){
    ND::THandle<ND::TTPCVolGroup> group1 = groups.at(i);
    for(int j=i+1; j<nGroups; ++j){
      ND::THandle<ND::TTPCVolGroup> group2 = groups.at(j);

      ND::THandle<ND::TTPCOrderedVolGroup> path = ConnectGroupPair(group1, group2);
      // make key, and vector of all redundancies from this pair
      std::pair<int, int> key = std::make_pair(i, j);
      std::vector<int> redundancies;

      // kill any groups associated with this path other than its direct ends
      for(int k=0; k<nGroups; ++k){
        if(k == i) continue;
        if(k == j) continue;
        ND::THandle<ND::TTPCVolGroup> group3 = groups.at(k);

        float minCost = 9999.;

        //bool conI = false;
        //bool conJ = false;
        for(std::vector<ND::TTPCPathVolume*>::iterator pathVolIt = path->begin(); pathVolIt != path->end(); ++pathVolIt){
          ND::TTPCPathVolume* pathVol = *pathVolIt;
          ND::TTPCUnitVolume* vol = pathVol->GetUnitVolume();

          float curCost = FindConnectionCost(group3, vol, true, true, true, edgeClearDist);
          minCost = std::min(minCost, curCost);
          /*if(curCost < edgeClearDist){
            // also check for CONSERVATIVE redundancies, to take precedence over non-conservative ones
            if(!conI){
              if(FindConnectionCost(group1, vol, true, true, true, edgeClearDist) >= edgeClearDist){
                conI = true;
              };
            };
            if(!conJ){
              if(FindConnectionCost(group2, vol, true, true, true, edgeClearDist) >= edgeClearDist){
                conJ = true;
              };
            };
          };*/
        };
        // if found, increment cost of group and add it to redundancy lists of connecting groups
        if(minCost < edgeClearDist){
          redundancies.push_back(k);

          /*if(conI && conJ){
            conservativeCrossRedundancies.at(i).push_back(k);
            conservativeCrossRedundancies.at(j).push_back(k);
          };*/
        };
      };
      // add list to map
      crossRedundancies[key] = redundancies;
    };
  };
  /*
  {
    std::cout << "\033[1;34m# " << "cross redundancies:  " << "\033[0;m" << std::endl;
    int i=0;
    for(std::map< std::pair<int, int>, std::vector<int> >::iterator circularEl = crossRedundancies.begin(); circularEl != crossRedundancies.end(); ++circularEl){
      std::vector<int> circular = circularEl->second;

      std::cout << "\033[1;35m# " << " " << circularEl->first.first << "-" << circularEl->first.second << ":  ( ";
      ++i;
      // remove all redundant elements not contained in the group
      for(std::vector<int>::iterator iIt = circular.begin(); iIt != circular.end(); ++iIt){
        std::cout << *iIt << " ";
      };
      std::cout << ")\033[0;m" << std::endl;
    };
  };*/

  // first add any and all circular redundancies; start by adding all groups
  std::vector< std::vector<int> > circularGroups;
  for(int i=0; i<nGroups; ++i){
    circularGroups.push_back( std::vector<int>(1, i) );
  };
  // now the fiddly bit - merge overlapping redundancies (ones on each others lists) until no more are to be merged
  bool maybeRedundant = true;
  while(maybeRedundant){
    maybeRedundant = false;

    for(std::vector< std::vector<int> >::iterator circular1It = circularGroups.begin(); circular1It != circularGroups.end(); ++circular1It){
      std::vector<int> circular1 = *circular1It;
      for(std::vector< std::vector<int> >::iterator circular2It = circular1It+1; circular2It != circularGroups.end(); ++circular2It){
        std::vector<int> circular2 = *circular2It;

        for(std::vector<int>::iterator iIt = circular1.begin(); iIt != circular1.end(); ++iIt){
          int i = *iIt;
          for(std::vector<int>::iterator jIt = circular2.begin(); jIt != circular2.end(); ++jIt){
            int j = *jIt;

            // vector of connections where i makees j redundant and vice versa
            std::vector<int> iInJTargets;
            std::vector<int> jInITargets;
            // look for i in any list spawned from j, or j in any from i
            for(std::map< std::pair<int, int>, std::vector<int> >::iterator crossRedundanciesEl = crossRedundancies.begin(); crossRedundanciesEl != crossRedundancies.end(); ++crossRedundanciesEl){
              std::pair<int, int> key = crossRedundanciesEl->first;
              std::vector<int> redundancies = crossRedundanciesEl->second;


              if((key.first == i) ^ (key.second == i)){
                if(std::find(redundancies.begin(), redundancies.end(), j) != redundancies.end()){
                  if(key.first == i) iInJTargets.push_back(key.second);
                  if(key.second == i) iInJTargets.push_back(key.first);
                };
              };
              if((key.first == j) ^ (key.second == j)){
                if(std::find(redundancies.begin(), redundancies.end(), i) != redundancies.end()){
                  if(key.first == j) jInITargets.push_back(key.second);
                  if(key.second == j) jInITargets.push_back(key.first);
                };
              };
            };
            bool redundant = false;
            // the two make each other redundant iff they make each other redundant with the same target
            for(std::vector<int>::iterator iInJTargetIt = iInJTargets.begin(); iInJTargetIt != iInJTargets.end(); ++iInJTargetIt)
            for(std::vector<int>::iterator jInITargetIt = jInITargets.begin(); jInITargetIt != jInITargets.end(); ++jInITargetIt){
              if(*iInJTargetIt == *jInITargetIt){
                redundant = true;
              };
            };

            // if they're on each others groups, merge the two
            if(redundant){
              for(std::vector<int>::iterator j2It = circular2.begin(); j2It != circular2.end(); ++j2It){
                int j2 = *j2It;

                if(std::find(circular1.begin(), circular1.end(), j2) == circular1.end()){
                  circular1It->push_back(j2);
                };
              };
              circularGroups.erase(circular2It);
              maybeRedundant = true;
              break;
            };
            if(maybeRedundant) break;
          };
          if(maybeRedundant) break;
        };
        if(maybeRedundant) break;
      };
      if(maybeRedundant) break;
    };
  };
  /*std::cout << "\033[1;34m# " << "circular groups:  " << "\033[0;m" << std::endl;
  for(std::vector< std::vector<int> >::iterator circularIt = circularGroups.begin(); circularIt != circularGroups.end(); ++circularIt){
    std::vector<int> circular = *circularIt;

    std::cout << "\033[1;35m# " << "  ( ";
    // remove all redundant elements not contained in the group
    for(std::vector<int>::iterator iIt = circular.begin(); iIt != circular.end(); ++iIt){
      std::cout << *iIt << " ";
    };
    std::cout << ")\033[0;m" << std::endl;
  };*/

  // now, make list of all free groups, initially containing everything
  std::vector<int> freeGroups;
  for(int i=0; i<nGroups; ++i){
    freeGroups.push_back(i);
  };
  // erase anything that's made redundant by elements of two circular groups not containing itself
  for(std::vector< std::vector<int> >::iterator circular1It = circularGroups.begin(); circular1It != circularGroups.end(); ++circular1It)
  for(std::vector< std::vector<int> >::iterator circular2It = circular1It; circular2It != circularGroups.end(); ++circular2It){
    std::vector<int> circular1 = *circular1It;
    std::vector<int> circular2 = *circular2It;

    for(std::vector<int>::iterator iIt = circular1.begin(); iIt != circular1.end(); ++iIt)
    for(std::vector<int>::iterator jIt = circular2.begin(); jIt != circular2.end(); ++jIt){
      int i = *iIt;
      int j = *jIt;
      std::pair<int, int> key;
      if(i < j){
        key = std::make_pair(i, j);
      }
      else if(i > j){
        key = std::make_pair(j, i);
      }
      else{
        continue;
      };
      std::vector<int> redundancies = crossRedundancies.at(key);

      for(std::vector<int>::iterator kIt = redundancies.begin(); kIt != redundancies.end(); ++kIt){
        int k = *kIt;

        // make sure it isn't in the circulars themselves
        if(std::find(circular1.begin(), circular1.end(), k) == circular1.end())
        if(std::find(circular2.begin(), circular2.end(), k) == circular2.end()){
          std::vector<int>::iterator findK = std::find(freeGroups.begin(), freeGroups.end(), k);
          if(findK != freeGroups.end()){
            // make sure it isn't in the two groups using for this check
            freeGroups.erase(findK);
          };
        };
      };
    };
  }
  // remove all elements not in the free list
  for(std::vector< std::vector<int> >::iterator circularIt = circularGroups.begin(); circularIt != circularGroups.end(); ++circularIt){
    std::vector<int> circular = *circularIt;
    std::vector<int>::iterator iIt = circular.begin();
    while(iIt != circular.end()){
      int i = *iIt;
      if(std::find(freeGroups.begin(), freeGroups.end(), i) == freeGroups.end()){
        iIt = circular.erase(iIt);
      }
      else{
        ++iIt;
      };
    };
  };
  // erase everything that's in a circular group and NOT the highest charge member
  for(std::vector< std::vector<int> >::iterator circularIt = circularGroups.begin(); circularIt != circularGroups.end(); ++circularIt){
    std::vector<int> circular = *circularIt;

    // remove all elements except one, by the following:
    //  a) <COMMENTED OUT> remove everything killed by the conservative redundancies list, unless that would be the whole group
    //  b) remove everything remaining but highest charge group

    // add everything in a non-shared conservative list
    /*
    std::vector<int> conservativeRedundancies;
    for(std::vector<int>::iterator iIt = circular.begin(); iIt != circular.end(); ++iIt)
    for(std::vector<int>::iterator jIt = iIt+1; jIt != circular.end(); ++jIt){
      int i = *iIt;
      int j = *jIt;

      std::vector<int> conI = conservativeCrossRedundancies.at(i);
      std::vector<int> conJ = conservativeCrossRedundancies.at(j);
      bool iInJ = (std::find(conI.begin(), conI.end(), j) != conI.end());
      bool jInI = (std::find(conJ.begin(), conJ.end(), i) != conJ.end());

      if(iInJ && !jInI){
        conservativeRedundancies.push_back(i);
      }
      else if(jInI && !iInJ){
        conservativeRedundancies.push_back(j);
      };
    };
    if(conservativeRedundancies.size() < circular.size()){
      // remove all elements found this way
      std::vector<int>::iterator kIt = circular.begin();
      while(kIt != circular.end()){
        int k = *kIt;
        std::vector<int>::iterator findK = std::find(freeGroups.begin(), freeGroups.end(), k);
        if(findK != freeGroups.end()){
          freeGroups.erase(findK);
        };
        if(std::find(conservativeRedundancies.begin(), conservativeRedundancies.end(), k) != conservativeRedundancies.end()){
          kIt = circular.erase(kIt);
        }
        else{
          ++kIt;
        };
      };
    };*/

    // now remove all but the highest charge element
    int leastRedundantID = -1;
    int leastRedundantCost = -9999.;
    for(std::vector<int>::iterator iIt = circular.begin(); iIt != circular.end(); ++iIt){
      int i = *iIt;
      float cost = groups.at(i)->GetAverageCharge();

      if(cost > leastRedundantCost){
        leastRedundantCost = cost;
        leastRedundantID = i;
      };
    };

    for(std::vector<int>::iterator iIt = circular.begin(); iIt != circular.end(); ++iIt){
      int i = *iIt;

      if(i != leastRedundantID){
        std::vector<int>::iterator findI = std::find(freeGroups.begin(), freeGroups.end(), i);
        if(findI != freeGroups.end()){
          freeGroups.erase(findI);
        };
      };
    };
  };
  // after all of this, if fewer than two ends remain try and add the most separated ones
  if(freeGroups.size() < 1){
    int maxI = -1;
    int maxJ = -1;
    float maxCost = -9999.;
    for(int i=0; i<nGroups; ++i){
      ND::THandle<ND::TTPCVolGroup> group1 = groups.at(i);
      for(int j=i+1; j<nGroups; ++j){
        ND::THandle<ND::TTPCVolGroup> group2 = groups.at(j);

        float cost = FindConnectionCost(group1, group2, true, true, true);
        if(cost > maxCost){
          maxCost = cost;
          maxI = i;
          maxJ = j;
        };
      };
    };
    if(maxI >= 0) freeGroups.push_back(maxI);
    if(maxJ >= 0) freeGroups.push_back(maxJ);
  }
  else if(freeGroups.size() < 2){
    int maxI = freeGroups.at(0);
    int maxJ = -1;
    float maxCost = -9999.;
    ND::THandle<ND::TTPCVolGroup> group1 = groups.at(maxI);
    for(int j=0; j<nGroups; ++j){
      if(j != maxI){
        ND::THandle<ND::TTPCVolGroup> group2 = groups.at(j);

        float cost = FindConnectionCost(group1, group2, true, true, true);
        if(cost > maxCost){
          maxCost = cost;
          maxJ = j;
        };
      };
    };
    if(maxJ >= 0) freeGroups.push_back(maxJ);
  };

  // and finally, add all free groups to out groups
  for(std::vector<int>::iterator iIt = freeGroups.begin(); iIt != freeGroups.end(); ++iIt){
    int i = *iIt;
    outGroups.push_back(groups.at(i));
  };

  return outGroups;
}
std::vector< ND::THandle<ND::TTPCVolGroup> > ND::TTPCAStar::ExperimentalClearVertexRedundancies(std::vector< ND::THandle<ND::TTPCVolGroup> > groups, int maxNo){
  // merge together vertices
  int vertexClearDist = fLayout->GetAltVertexHitConnectDist();
  std::vector< ND::THandle<ND::TTPCVolGroup> > mergedGroups = MergeGroupsAStar(groups, vertexClearDist);

  return mergedGroups;
}

std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > ND::TTPCAStar::ClearVertexConnectionRedundancies(ND::TTPCVolGroupMan* volGroupMan, std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > paths, std::vector< ND::THandle<ND::TTPCVolGroup> > vertices){
  std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > survivors;

  // start with all groups as survivors
  for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator path = paths.begin(); path != paths.end(); ++path) survivors.push_back(*path);

  for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator path = paths.begin(); path != paths.end(); ++path){
    for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator vert = vertices.begin(); vert != vertices.end(); ++vert){
      bool redundancy = volGroupMan->GetPathVolOverlap(*path, (*vert)->GetAverageVol(), ND::TTPCConnection::vertexPath);
      if(redundancy){
        for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator it = survivors.begin(); it != survivors.end(); ++it)
          if(*it == *path){
            survivors.erase(it);
            break;
          };
        break;
      };
    };
  };

  return survivors;
}
ND::THandle<ND::TTPCOrderedVolGroup> ND::TTPCAStar::ConnectGroupPair(ND::THandle<ND::TTPCVolGroup> group1, ND::THandle<ND::TTPCVolGroup> group2, bool vertexGroup1, bool vertexGroup2, bool fullASICPenalty, bool extendMode){
  ND::THandle<ND::TTPCOrderedVolGroup> connection(new ND::TTPCOrderedVolGroup(fLayout));
  // get unique ids and A* indices
  ND::TTPCUnitVolume* startCell = group1->GetAverageVol();
  ND::TTPCUnitVolume* endCell = group2->GetAverageVol();

  if(!startCell || !endCell) return connection;

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
    return connection;
  };

  // add front and back groups to the returned group (front corresponds to hits near size()-1 index, back to hits near 0 index)
  connection->AddFrontHits(group1);
  connection->AddBackHits(group2);
  // set status of front and back groups as vertices (or not)
  connection->SetFrontIsVertex(vertexGroup1);
  connection->SetBackIsVertex(vertexGroup2);

  // start at end index and add the chain of cells connecting it from the start to the connection
  ND::TTPCAStarPoint* curPoint = endPoint;
  int nPoints = fAStarPoints.size();
  for(int i=0; i < nPoints; i++){
    connection->AddCell(curPoint->vol);

    ND::TTPCAStarPoint* nextPoint = curPoint->aStarParent;
    if(!nextPoint) break; 
    curPoint = nextPoint;
  };
  // add the first cell so long as the connection isn't one cell long
  if(curPoint != endPoint) connection->AddCell(curPoint->vol);

  return connection;
}
float ND::TTPCAStar::FindConnectionCost(ND::THandle<ND::TTPCVolGroup> group1, ND::THandle<ND::TTPCVolGroup> group2, bool fullASICPenalty, bool extendMode, bool reduced, float maxCost){
  ND::TTPCUnitVolume* vol1 = group1->GetAverageVol();
  ND::TTPCUnitVolume* vol2 = group2->GetAverageVol();

  return FindConnectionCost(vol1, vol2, fullASICPenalty, extendMode, reduced, maxCost);
}
float ND::TTPCAStar::FindConnectionCost(ND::THandle<ND::TTPCVolGroup> group, ND::TTPCUnitVolume* vol, bool fullASICPenalty, bool extendMode, bool reduced, float maxCost){
  ND::TTPCUnitVolume* groupVol = group->GetAverageVol();

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

void ND::TTPCAStar::AssociateBestHits(ND::TTPCVolGroupMan* volGroupMan, std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > inPaths, bool fullASICPenalty, bool extendMode, float maxCost){
  // add hits in path to group, making sure the same isn't added twice
  std::vector< ND::THandle<ND::TTPCVolGroup> > hitsGroups;
  std::set<ND::TTPCUnitVolume*> volsAdded;
  for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator inPathIt = inPaths.begin(); inPathIt != inPaths.end(); ++inPathIt){
    ND::THandle<ND::TTPCOrderedVolGroup> inPath = *inPathIt;

    ND::THandle<ND::TTPCVolGroup> hitsGroup (new ND::TTPCVolGroup(fLayout));
    for(std::vector<ND::TTPCPathVolume*>::iterator pathVolIt = inPath->begin(); pathVolIt != inPath->end(); ++pathVolIt){
      ND::TTPCPathVolume* pathVol = *pathVolIt;
      ND::TTPCUnitVolume* vol = pathVol->GetUnitVolume();

      if(!volsAdded.count(vol)){
        hitsGroup->AddHit(vol);
        volsAdded.insert(vol);
      };
    };

    hitsGroups.push_back(hitsGroup);
  };

  // merge best hits in groups
  MergeBestHits(volGroupMan, hitsGroups, fullASICPenalty, extendMode, maxCost);

  // add merged hits to appropriate paths
  for(unsigned int i=0; i<inPaths.size(); ++i){
    ND::THandle<ND::TTPCOrderedVolGroup> inPath = inPaths.at(i);
    ND::THandle<ND::TTPCVolGroup> hitsGroup = hitsGroups.at(i);

    inPath->AddExtendedHits(hitsGroup);
  };
}
void ND::TTPCAStar::MergeBestHits(ND::TTPCVolGroupMan* volGroupMan, ND::THandle<ND::TTPCVolGroup> inGroup, bool fullASICPenalty, bool extendMode, float maxCost){
  std::vector< ND::THandle<ND::TTPCVolGroup> > dummyVector;
  dummyVector.push_back(inGroup);

  MergeBestHits(volGroupMan, dummyVector, fullASICPenalty, extendMode, maxCost);
}
void ND::TTPCAStar::MergeBestHits(ND::TTPCVolGroupMan* volGroupMan, std::vector< ND::THandle<ND::TTPCVolGroup> > inGroups, bool fullASICPenalty, bool extendMode, float maxCost){
  // reset hits
  RebootHits();

  // define container of all near hits currently under consideration, starting with path hits
  std::vector<ND::TTPCAStarPoint*> openSet;
  std::vector<ND::TTPCAStarPoint*> definingSet (inGroups.size(), 0);
  int nGroup = 0;
  for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator groupIt = inGroups.begin(); groupIt != inGroups.end(); ++groupIt){
    ND::THandle<ND::TTPCVolGroup> group = *groupIt;

    for(std::map<long, ND::TTPCUnitVolume*>::iterator hitEl = group->begin(); hitEl != group->end(); ++hitEl){
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
  for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator groupIt = inGroups.begin(); groupIt != inGroups.end(); ++groupIt){
    ND::THandle<ND::TTPCVolGroup> group = *groupIt;
    ND::TTPCAStarPoint* definingPoint = definingSet[nGroup];
    ++nGroup;

    if(definingPoint != 0){
      for(std::vector<ND::TTPCAStarPoint*>::iterator pntIt = fAStarPoints.begin(); pntIt != fAStarPoints.end(); ++pntIt){
        ND::TTPCAStarPoint* pnt = *pntIt;

        // check if associated with this path
        if(pnt->aStarParent == definingPoint){
          // add vol to group
          group->AddHit(pnt->vol);
        };
      };
    };
  };
}

void ND::TTPCAStar::GetNearHitConnections(ND::TTPCVolGroupMan* volGroupMan, ND::TTPCAStarPoint* point, bool extendMode){
  // only search for friends if none are already found
  if(!point->aStarHasFriends){
    // volume group to check for nearby hits
    ND::THandle<ND::TTPCVolGroup> volCheck(new ND::TTPCVolGroup(fLayout));
    volCheck->AddHitMap(fHitMap);

    // find near hits and (friends)
    ND::THandle<ND::TTPCVolGroup> nearHits = volGroupMan->GetNearHits(volCheck, point->vol, ND::TTPCConnection::path);

    // save cost of connecting to each of those friends
    for(std::map<long, ND::TTPCUnitVolume*>::iterator it = nearHits->begin(); it != nearHits->end(); ++it){
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
    ND280Warn("Path not found between two points");
    return 1;
  };
  return 0;
}
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
}

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
