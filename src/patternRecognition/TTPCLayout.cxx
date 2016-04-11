// eddy
#include "TTPCLayout.hxx"

trex::TTPCLayout::TTPCLayout(){

  //MDH
  //This all needs reimplementing, and we need to think about if we really need
  //all these parameters...
  /*
  // values for cuts on charge
  fChargeCut = trex::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.PatRec.ChargeCut");
  fEarlyNegativePeakCut = trex::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.PatRec.EarlyNegativePeakCut");
  fLateNegativePeakCut = trex::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.PatRec.LateNegativePeakCut");
  fASICSaturationCut = trex::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.PatRec.ASICSaturationCut");
  fASICSubOccupancyCut = trex::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.PatRec.ASICSubOccupancyCut");
  fASICOccupancyCut = trex::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.PatRec.ASICOccupancyCut");
  fASICSatExpansion = trex::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.PatRec.ASICSatExpansion");
  fASICOccExpansion = trex::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.PatRec.ASICOccExpansion");
  fASICSplittingY = trex::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.PatRec.ASICSplittingY");
  fASICSplittingZ = trex::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.PatRec.ASICSplittingZ");

  // sampling time and drift speeds
  fDriftSpeed = trex::tpcCalibration().GetDriftVelocity();

  // set pre-defined geometry variables
  fPadGap = trex::TGeomInfo::Get().TPC().GetPadGap();
  fPadPitchY = trex::TGeomInfo::Get().TPC().GetPadYPitch();
  fPadPitchZ = trex::TGeomInfo::Get().TPC().GetPadXPitch();

  fMMPads = trex::TGeomInfo::Get().TPC().GetPadCount();
  fYPads = trex::TGeomInfo::Get().TPC().GetPadRows();
  fZPads = trex::TGeomInfo::Get().TPC().GetPadColumns();

  // set x cell geometry variables
  fXCellSize = trex::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.PatRec.XCellSize") * fPadPitchZ; 
  fTWidth = fXCellSize / fDriftSpeed;

  // set if groups are broken in x, y and z
  fJumpX = bool( trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.JumpX") );
  fJumpY = bool( trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.JumpY") );
  fJumpZ = bool( trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.JumpZ") );

  // set offsets
  fGapOffsetX = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.GapOffsetX");
  fGapOffsetY = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.GapOffsetY");
  fGapOffsetZ = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.GapOffsetZ");
  fGapOffsetAdjacent = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.GapOffsetAdjacent");

  // set gap variables
  fXShiftFromC = fGapOffsetX;
  fYShiftFromXZ = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.YShiftFromXZ");
  fYShiftFromMM = fYPads + fGapOffsetY;
  fZShiftFromMM = fZPads + fGapOffsetZ;
  fZShiftFromTPC = (fZPads + fGapOffsetZ)*3;

  // set useable pattern and path sizes
  fMinPatternPads = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.MinPatternPads");
  fMinPathClusters = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.MinPathClusters");

  // set variables for preliminary delta ray search
  fDeltaSpreadRate = trex::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.PatRec.DeltaSpreadRate");

  // set variables for primary track edge search
  fUseIndirectEdges = (bool)trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.UseIndirectEdges");
  fEdgeLayers = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.EdgeLayers");

  // set pattern recognition and path finding connection distances
  fConnectDistX = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.ConnectDistX");
  fConnectDistY = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.ConnectDistY");
  fConnectDistZ = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.ConnectDistZ");

  // set charactaristic structure distances
  fStructDistX = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.StructDistX");
  fStructDistY = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.StructDistY");
  fStructDistZ = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.StructDistZ");

  // characterising distances based on HitDist
  fAStarConnectDist = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.AStarConnectDist");
  fPathHitConnectDist = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.PathHitConnectDist");
  fExtraHitConnectDist = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.ExtraHitConnectDist");

  // characteristic distances based on StructDist
  fEdgeMergeStructDist = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.EdgeMergeStructDist");
  fVertexMergeStructDist = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.VertexMergeStructDist");
  fVertexFindStructDist = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.VertexFindStructDist");
  fVertexHitStructDist = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.VertexHitStructDist");
  fVertexPathStructDist = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.VertexPathStructDist");
  fClusterMergeStructDist = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.ClusterMergeStructDist");

  // variables used when fUseAltHitAssociation is true
  fAltPathHitConnectDist = trex::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.PatRec.AltPathHitConnectDist");
  fAltExtraHitConnectDist = trex::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.PatRec.AltExtraHitConnectDist");
  fAltEdgeHitConnectDist = trex::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.PatRec.AltEdgeHitConnectDist");
  fAltVertexHitConnectDist = trex::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.PatRec.AltVertexHitConnectDist");

  // variables used to check whether to count something as an x path
  fXPathMaxPads = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.XPathMaxPads");
  fXPathMinEndRatio = trex::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.PatRec.XPathMinEndRatio");

  // set edge parameters
  fEdgePreDist = trex::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.PatRec.EdgePreDist");
  fEdgePreAng = trex::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.PatRec.EdgePreAng");
  fEdgeMinDist = trex::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.PatRec.EdgeMinDist");
  fEdgeMinDistLow = trex::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.PatRec.EdgeMinDistLow");
  fEdgeMinHits = trex::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.PatRec.EdgeMinHits");
  fEdgeMaxSigma = trex::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.PatRec.EdgeMaxSigma");
  fEdgeRange = trex::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.PatRec.EdgeRange");
  fEdgeThreshold = trex::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.PatRec.EdgeThreshold");
  fEdgeOffset = trex::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.PatRec.EdgeOffset");

  // set parameters for forming HV clusters
  fHVClusterExtrapolateDist = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.HVClusterExtrapolateDist");
  fHVClusterExtrapolateLimit = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.HVClusterExtrapolateLimit");
  fClusterConnectDist = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.ClusterConnectDist");
  fMergeDist = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.MergeDist");
  fHVEdgeDist = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.HVEdgeDist");
  fThresholdAngleRange = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.ThresholdAngleRange");
  fThresholdAngle = trex::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.PatRec.ThresholdAngle");
  fDichotomyCutoff = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.DichotomyCutoff");

  // set parameters relating to features along paths
  fHVClusterMaxIso = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.HVClusterMaxIso");

  // set parameters for breaking long x clusters
  fXSizeThreshold = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.XSizeThreshold");
  fPathSizeThreshold = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.PathSizeThreshold");
  fBreakInMiddle = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.BreakInMiddle");

  // set parameters relating to cleaning up hits near junctions
  fAnomCheckDist = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.AnomCheckDist");
  fAnomProjectDist = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.AnomProjectDist");
  fAnomMaxOffs = trex::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.PatRec.AnomMaxOffs");

  // set parameters relating to cleaning up hits near junctions
  fAnomCheckDist = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.AnomCheckDist");
  fAnomProjectDist = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.AnomProjectDist");
  fAnomMaxOffs = trex::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.PatRec.AnomMaxOffs");

  // defualt scales and heuristic factor for A* algorithm
  fAStarXScale = trex::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.PatRec.AStarXScale");
  fAStarYScale = trex::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.PatRec.AStarYScale");
  fAStarZScale = trex::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.PatRec.AStarZScale");

  fAStarHeuristicFactor = trex::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.PatRec.AStarHeuristicFactor");
  fAStarPathologyPenalty = trex::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.PatRec.AStarPathologyPenalty");
  fAStarAssociatePathologyPenalty = trex::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.PatRec.AStarAssociatePathologyPenalty");

  // minimum fractional and absolute non-delta hits for avoiding EM classification
  fNonDeltaFraction= trex::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.PatRec.NonDeltaFraction");
  fNonDelta = trex::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatRec.NonDelta");
  
  // other default values
  fTNMin = 0.;
  fTNMax = -1.;
  fTPMin = 0.;
  fTPMax = -1.;
  fTNegativeBins = 0;
  fTPositiveBins = 0;
  fTBins = 0;*/
}
trex::TTPCLayout::~TTPCLayout(){
}

void trex::TTPCLayout::SetTimeRanges(double tNMin, double tNMax, double tPMin, double tPMax){
  fTNMin = tNMin;
  fTNMax = tNMax;
  fTPMin = tPMin;
  fTPMax = tPMax;

  bool binsN = tNMin < tNMax; // do hits exist in negative half?
  bool binsP = tPMin < tPMax; // do hits exist in positibe half?
  fXCathodeCross = (binsN && binsP);

  if(binsN){
    fTNegativeBins = (int)( (fTNMax-fTNMin) / fTWidth );
  }
  else{
    fTNegativeBins = 0;
  };
  if(binsP){
    fTPositiveBins = (int)( (fTPMax-fTPMin) / fTWidth );
  }
  else{
    fTPositiveBins = 0;
  };
  fTBins = fTNegativeBins + fTPositiveBins;
  if(fXCathodeCross) fTBins += fXShiftFromC;
}
void trex::TTPCLayout::SetRanges(int minX,int maxX, int minY,int maxY, int minZ,int maxZ){
  // set x, y and z minima, maxima and sizes
  fMinX = minX;
  fMinY = minY;
  fMinZ = minZ;

  fMaxX = maxX;
  fMaxY = maxY;
  fMaxZ = maxZ;

  // total number of cells in x, y and z is max - min, plus one extra cell to be inclusive
  fSizeX = fMaxX - fMinX + 1;
  fSizeY = fMaxY - fMinY + 1;
  fSizeZ = fMaxZ - fMinZ + 1;
}
void trex::TTPCLayout::GetRanges(int& sizeX,int& minX,int& maxX, int& sizeY,int& minY,int& maxY, int& sizeZ,int& minZ,int& maxZ){
  // get x, y and z minima, maxima and sizes
  sizeX = fSizeX;
  minX = fMinX;
  maxX = fMaxX;

  sizeY = fSizeY;
  minY = fMinY;
  maxY = fMaxY;

  sizeZ = fSizeZ;
  minZ = fMinZ;
  maxZ = fMaxZ;
}
void trex::TTPCLayout::GetRanges(int& sizeX,int& minX,int& maxX, int& sizeY,int& minY,int& maxY, int axis){
  // get 2D minima, maxima and sizes in x, y or z view (with axis of 1, 2 or 3 respectively)
  switch(axis){
    case 1:
      sizeX = fSizeZ;
      minX = fMinZ;
      maxX = fMaxZ;

      sizeY = fSizeY;
      minY = fMinY;
      maxY = fMaxY;

      break;
    case 2:
      sizeX = fSizeZ;
      minX = fMinZ;
      maxX = fMaxZ;

      sizeY = fSizeX;
      minY = fMinX;
      maxY = fMaxX;

      break;
    case 3:
      sizeX = fSizeY;
      minX = fMinY;
      maxX = fMaxY;

      sizeY = fSizeX;
      minY = fMinX;
      maxY = fMaxX;

      break;
  };
}

long trex::TTPCLayout::Mash(int x, int y, int z){
  // convert x, y and z id to unique id
  long valX = long(x - fMinX);
  long valY = long(y - fMinY);
  long valZ = long(z - fMinZ);
  return (valZ * fSizeX * fSizeY) + (valY * fSizeX) + valX;
}
long trex::TTPCLayout::MashYZ(int y, int z){
  // convert x, y and z id to unique id
  long valY = long(y - fMinY);
  long valZ = long(z - fMinZ);
  return (valZ * fSizeY) + valY;
}
long trex::TTPCLayout::SafeMash(int x, int y, int z){
  // return -1 if x, y or z id are invalid
  if (x < fMinX || x > fMaxX) return -1;
  if (y < fMinY || y > fMaxY) return -1;
  if (z < fMinZ || z > fMaxZ) return -1;
  // convert x, y and z id to unique id
  return Mash(x, y, z);
}
trex::TTPCCell3D trex::TTPCLayout::UnMash(long id){
  // convert unique id to x, y and z id
  trex::TTPCCell3D cell;

  cell.x = id % fSizeX;
  id = (id - cell.x) / fSizeX;
  cell.y = id % fSizeY;
  id = (id - cell.y) / fSizeY;
  cell.z = id;

  cell.x += fMinX;
  cell.y += fMinY;
  cell.z += fMinZ;

  return cell;
}

trex::TTPCPadStruct trex::TTPCLayout::GlobalXYZToPos(TVector3 pos){
  trex::TTPCPadStruct padStruct;

  //MDH
  //This needs reimplementing with new geometry lookup
  //or just simplifying
  // get pad geometry information
  /*
  trex::TGeometryId id;
  trex::TGeomInfo::Get().TPC().GlobalXYZToGeomId(pos, id);

  // extract location for tpc, gas half, mm unit and pad
  padStruct.tpc = trex::TGeomInfo::Get().TPC().GeomIdToTPC(id);
  padStruct.half = trex::TGeomInfo::Get().TPC().GeomIdToHalf(id);
  padStruct.mm = trex::TGeomInfo::Get().TPC().GeomIdToMM(id);
  padStruct.rawpad = trex::TGeomInfo::Get().TPC().GeomIdToPad(id);

  // shift to consistent pattern of MM pads (consistent order, consistent ids in x, y and z)
  padStruct.pad = padStruct.rawpad;
  if (padStruct.mm < 6){
    padStruct.pad = (fMMPads-1) - padStruct.pad;
  };
  if(padStruct.half == 0){
    int cpad = padStruct.pad%fYPads;
    padStruct.pad -= cpad;
    padStruct.pad = (fMMPads-fYPads) - padStruct.pad;
    padStruct.pad += cpad;
  };

  // shift to consistent pattern of MM module ids
  if(padStruct.half == 1){
    padStruct.mm = (padStruct.mm + 6) % 12;
    };*/

  // mm numbers and mm pad numbers should now be oriented consistently across all modules
  return padStruct;
}

//MDH
//Got rid of this sig since it uses a THit. Replace invocation of this sig in TTPCTRExPatAlgorithm
//with the one below

/*
trex::TTPCCellInfo3D trex::TTPCLayout::GetPadPosID(trex::THandle<trex::THit> hit, int tpcMask){
  double time;
  trex::THandle<trex::TTPCHitPad> hitPad = hit;

  if(hitPad){
    std::vector<double> peakTimes = hitPad->GetPeakTimes();
    if(peakTimes.size()){
      time = 0.;
      for(std::vector<double>::iterator peakTimeIt = peakTimes.begin(); peakTimeIt != peakTimes.end(); ++peakTimeIt){
        time += *peakTimeIt;
      };
      time /= (double)peakTimes.size();
    }
    else{
      time = hitPad->GetTime();
    };
  }
  else{
    time = hit->GetTime();
  };

  return GetPadPosID(hit->GetPosition(), time, tpcMask);
}
*/

trex::TTPCCellInfo3D trex::TTPCLayout::GetPadPosID(TVector3 pos, double time, int tpcMask){
  trex::TTPCCellInfo3D cellInfo;

  // get pad id
  trex::TTPCPadStruct pad = GlobalXYZToPos(pos);

  // initialise edge status to 0
  cellInfo.edgeX = 0;
  cellInfo.edgeY = 0;
  cellInfo.edgeZ = 0;
  if (tpcMask > 0){
    // if tpc mask is on and tpc is not that specified, pad is invalid
    if (pad.tpc+1 != tpcMask) pad.tpc = -1;
  };
  if ((pad.tpc < 0) || (pad.half < 0) || (pad.mm < 0) || (pad.pad < 0)){
    // if pad is invalid, set all ids to -1 and return
    cellInfo.x = -1;
    cellInfo.y = -1;
    cellInfo.z = -1;
    return cellInfo;
  };

  // calculate cell x by finding a time bin depending on the cathode half and the hit time
  if(pos.X() < 0){
    // negative side - just bin normally 
    cellInfo.x = (int)( (time-fTNMin) / fTWidth );
  }
  else{
    // positive side - subtract bin number from maximum bin
    cellInfo.x = fTBins - (int)( (time-fTPMin) / fTWidth );
  };

  // get cell raw y
  cellInfo.y = pad.pad%fYPads;
  // get cell raw z
  cellInfo.z = (pad.pad-cellInfo.y)/fYPads;

  // if next to central cathode, set x edges appropriately 
  if (fXCathodeCross){
    if(pos.X() < 0){
      if( (fTNMax-time) < fTWidth ) cellInfo.edgeX = 1;
    }
    else{
      if( (fTPMax-time) < fTWidth ) cellInfo.edgeX = -1;
    };
  };

  // if at y mm gap, set y edge appropriately
  if (cellInfo.y == fYPads - 1) cellInfo.edgeY = 1;
  else if(cellInfo.y == 0) cellInfo.edgeY = -1;

  // if at z mm gap, set z edge appropriately
  if (cellInfo.z == fZPads - 1) cellInfo.edgeZ = 1;
  else if(cellInfo.z == 0) cellInfo.edgeZ = -1;

  // if pad is dead by default, set edges to -2
  if (pad.rawpad == 0 || pad.rawpad == 1){
    cellInfo.edgeX = -2;
    cellInfo.edgeY = -2;
    cellInfo.edgeZ = -2;
  };

  // shift y id down if in a shifted mm column
  if (pad.mm < 6){
    if (pad.half == 0) cellInfo.y += fYShiftFromXZ;
  }
  else{
    if (pad.half == 1) cellInfo.y += fYShiftFromXZ;
    cellInfo.z += fZShiftFromMM;
  };
  // increase y by position of mm volume in row times number of pads contained
  cellInfo.y += (5 - (pad.mm%6)) * fYShiftFromMM;
  // increase z by position of tpc times number of pads in tpc
  cellInfo.z += pad.tpc * fZShiftFromTPC;

  // set segments appropriately depending on MM position
  cellInfo.segX = pad.half;
  cellInfo.segY = pad.mm % 6;
  cellInfo.segZ = 2*pad.tpc + ( (int)pad.mm < 6 );

  //if (cellInfo.z > 1000)std::cout << "\033[1;33m# " << "z: " << cellInfo.z << "\033[0;m" << std::endl;
  return cellInfo;
}

void trex::TTPCLayout::GetTypeDistances(int& distX, int& distY, int& distZ, trex::TTPCConnection::Type type){
  if (type==trex::TTPCConnection::path){
    distX = fConnectDistX * fAStarConnectDist;
    distY = fConnectDistY * fAStarConnectDist;
    distZ = fConnectDistZ * fAStarConnectDist;
  }
  else if (type==trex::TTPCConnection::pathHits){
    distX = fConnectDistX * fPathHitConnectDist;
    distY = fConnectDistY * fPathHitConnectDist;
    distZ = fConnectDistZ * fPathHitConnectDist;
  }
  else if (type==trex::TTPCConnection::extraHits){
    distX = fConnectDistX * fExtraHitConnectDist;
    distY = fConnectDistY * fExtraHitConnectDist;
    distZ = fConnectDistZ * fExtraHitConnectDist;
  }
  else if (type==trex::TTPCConnection::edgeMerge){
    distX = fStructDistX * fEdgeMergeStructDist;
    distY = fStructDistY * fEdgeMergeStructDist;
    distZ = fStructDistZ * fEdgeMergeStructDist;
  }
  else if (type==trex::TTPCConnection::vertexMerge){
    distX = fStructDistX * fVertexMergeStructDist;
    distY = fStructDistY * fVertexMergeStructDist;
    distZ = fStructDistZ * fVertexMergeStructDist;
  }
  else if (type==trex::TTPCConnection::vertexFind){
    distX = fStructDistX * fVertexFindStructDist;
    distY = fStructDistY * fVertexFindStructDist;
    distZ = fStructDistZ * fVertexFindStructDist;
  }
  else if (type==trex::TTPCConnection::vertexHit){
    distX = fStructDistX * fVertexHitStructDist;
    distY = fStructDistY * fVertexHitStructDist;
    distZ = fStructDistZ * fVertexHitStructDist;
  }
  else if (type==trex::TTPCConnection::vertexPath){
    distX = fStructDistX * fVertexPathStructDist;
    distY = fStructDistY * fVertexPathStructDist;
    distZ = fStructDistZ * fVertexPathStructDist;
  }
  else if (type==trex::TTPCConnection::clusterMerge){
    distX = fStructDistX * fClusterMergeStructDist;
    distY = fStructDistY * fClusterMergeStructDist;
    distZ = fStructDistZ * fClusterMergeStructDist;
  };
}
