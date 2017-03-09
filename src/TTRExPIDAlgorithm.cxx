#include "TTRExPIDAlgorithm.hxx" 

//root
#include <TVector3.h>

//********************************************************************
trex::TTRExPIDAlgorithm::TTRExPIDAlgorithm( ){};

//trex::TTRExPIDAlgorithm::~TTRExPIDAlgorithm( ){};
//********************************************************************

void trex::TTRExPIDAlgorithm::Process(std::vector<trex::TTRExPattern>& allPatterns){
  
  for (auto patit=allPatterns.begin(); patit != allPatterns.end(); patit++) {
    trex::TTRExPattern& Pattern=*patit;
    std::vector<trex::TTRExPath>& Paths = Pattern.GetPaths();
    
    for(auto pathit=Paths.begin(); pathit !=Paths.end(); pathit++){
      trex::TTRExPath& path=*pathit;
      
      dEdx(path);      
    
    }
  }
}

void trex::TTRExPIDAlgorithm::dEdx(trex::TTRExPath& path){

  ChargeSum(path); 
  TrackLength(path);

  double charge = path.GetChargeSum();
  double length = path.GetTrackLength();

  double dEdx = charge/length;

  std::cout << "Calculated dE/dx as: " << dEdx << std::endl;

  path.SetdEdx(dEdx);
}



void trex::TTRExPIDAlgorithm::ChargeSum(trex::TTRExPath& path){

  double charge = 0;

  std::vector<trex::TTRExHVCluster*>& Clusters = path.GetClusters();
  
  for(auto clusterit=Clusters.begin(); clusterit != Clusters.end(); clusterit++){
    
    trex::TTRExHVCluster& Cluster=**clusterit;
    
    std::vector<trex::TTPCHitPad*>& Hits = Cluster.GetClusterHits();

    for(auto hitit=Hits.begin(); hitit != Hits.end(); hitit++){
    
      trex::TTPCHitPad& Hit = **hitit;
      
      charge += Hit.GetCharge();
    }
  }

  path.SetChargeSum(charge);

}


void trex::TTRExPIDAlgorithm::TrackLength(trex::TTRExPath& path){

  std::vector<trex::TTRExHVCluster*>::iterator tmpClu = path.GetClusters().begin();
  trex::TTRExHVCluster* FirstCluster = *tmpClu;
  
  TVector3 trackStart = FirstCluster->GetPosition();

  std::vector<trex::TTRExHVCluster*>::reverse_iterator tmpClu2 = path.GetClusters().rbegin();
  trex::TTRExHVCluster* LastCluster = *tmpClu2;

  TVector3 trackEnd = LastCluster->GetPosition();

  TVector3 track = trackEnd-trackStart;
  double trackLength = track.Mag();

}
