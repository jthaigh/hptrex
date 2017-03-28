// eddy
#include "TTPCLayout.hxx"

trex::TTPCLayout::TTPCLayout(){

  //MDH
  //This all needs reimplementing, and we need to think about if we really need
  //all these parameters...

  fHaveX=true;

  fBField=1.;

  //MDH
  //PLACEHOLDER!!!
  // sampling time and drift speeds
  fPadPitchY = 1.;
  fPadPitchZ = 1.;
  fXCellSize = 1. * fPadPitchZ; 

  // set useable pattern and path sizes
  fMinPatternPads = 4;
  fMinPathClusters = 4;

  fEdgeLayers = 2;

  // set pattern recognition and path finding connection distances
  //P.D. 
  //fConnectDistX = 3;
  fConnectDistX = 2;
  fConnectDistY = 2;
  fConnectDistZ = 2;

  // set charactaristic structure distances
  fStructDistX = 1;
  fStructDistY = 1;
  fStructDistZ = 1;

  // characterising distances based on HitDist
  fAStarConnectDist = 2;
  fPathHitConnectDist = 10;
  fExtraHitConnectDist = 4;

  // characteristic distances based on StructDist
  fEdgeMergeStructDist = 3;
  fVertexMergeStructDist = 6;
  fVertexFindStructDist = 3;
  fVertexHitStructDist = 2;
  fVertexPathStructDist = 6;
  fClusterMergeStructDist = 2;

  // variables used to check whether to count something as an x path
  fXPathMaxPads = 5;
  fXPathMinEndRatio = 2.;

  // set edge parameters
  fEdgePreDist = 8.;
  fEdgePreAng = 80.;
  fEdgeMinDist = 8.;
  fEdgeMinDistLow = 2.;
  fEdgeMinHits = 0.2;
  fEdgeMaxSigma = 0.3;
  fEdgeRange = 0.56;
  fEdgeThreshold = 0.52;
  fEdgeOffset = 0.3;

  // set parameters for forming HV clusters
  fHVClusterExtrapolateDist = 3;
  fHVClusterExtrapolateLimit = 20;
  fClusterConnectDist = 5;
  fMergeDist = 8;
  fHVEdgeDist = 12;
  fThresholdAngleRange = 3;
  fThresholdAngle = 55.;
  fDichotomyCutoff = 8;

  // set parameters relating to features along paths
  fHVClusterMaxIso = 3;

  // set parameters for breaking long x clusters
  fXSizeThreshold = 4;
  fPathSizeThreshold = 3;
  fBreakInMiddle = 0;

  // set parameters relating to cleaning up hits near junctions
  fAnomCheckDist = 2;
  fAnomProjectDist = 8;
  fAnomMaxOffs = 4.;

  // defualt scales and heuristic factor for A* algorithm
  fAStarXScale = 1.;
  fAStarYScale = 1.;
  fAStarZScale = 1.;

  fAStarHeuristicFactor = 1.2;
  
}
trex::TTPCLayout::~TTPCLayout(){
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

long trex::TTPCLayout::Mash(int x, int y, int z){
  // convert x, y and z id to unique id
  long valX = long(x - fMinX);
  long valY = long(y - fMinY);
  long valZ = long(z - fMinZ);
  return (valZ * fSizeX * fSizeY) + (valY * fSizeX) + valX;
}

long trex::TTPCLayout::SafeMash(int x, int y, int z){
  // return -1 if x, y or z id are invalid
  if (x < fMinX || x > fMaxX) return -1;
  if (y < fMinY || y > fMaxY) return -1;
  if (z < fMinZ || z > fMaxZ) return -1;
  // convert x, y and z id to unique id
  return Mash(x, y, z);
}

trex::TTPCCellInfo3D trex::TTPCLayout::GetPadPosID(TVector3 pos, int tpcMask){
  trex::TTPCCellInfo3D cellInfo;

  //MDH
  //PLACEHOLDER!!!
  cellInfo.x=pos.X();
  cellInfo.y=pos.Y();
  cellInfo.z=pos.Z();

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

