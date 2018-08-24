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
      TrackCleanliness(path);
      TrueTrackLength(path);
      PID(path);
    }
  }
}

void trex::TTRExPIDAlgorithm::dEdx(trex::TTRExPath& path){

  ChargeSum(path); 
  TrackLength(path);

  double charge = path.GetChargeSum();
  double length = path.GetTrackLength();

  double dedx;
  if(length!=0){
    dedx = charge/length;
  }else{dedx = charge;}

  std::cout << "Calculated dE/dx as: " << dedx << std::endl;

  path.SetdEdx(dedx);

  std::cout << "Path dE/dx was set" << std::endl;
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
  //std::cout << "CLUSTER POSITION 1: " << trackStart.Print() << std::endl;


  std::vector<trex::TTRExHVCluster*>::reverse_iterator tmpClu2 = path.GetClusters().rbegin();
  trex::TTRExHVCluster* LastCluster = *tmpClu2;

  TVector3 trackEnd = LastCluster->GetPosition();
  //std::cout <<"CLUSTER POSITION 2: " << trackEnd.Print() << std::endl;


  TVector3 track = trackEnd-trackStart;
  double mag = track.Mag();
  double trackLength = mag*2.34;
  
  std::cout << "*************************************************" << std::endl;
  std::cout << "TRACK LENGTH WAS: " << trackLength << std::endl; 
  std::cout << "*************************************************" << std::endl;

  path.SetTrackLength(trackLength);

}


//implement
void trex::TTRExPIDAlgorithm::TrueTrackLength(trex::TTRExPath& path){

  TVector3 initialPos =  path.GetInitialPosition();
  TVector3 finalPos = path.GetFinalPosition();

  TVector3 trueTrack = finalPos-initialPos;
  double trueLength = trueTrack.Mag();

  std::cout << "*************************************************" << std::endl;
  std::cout << "TRUE TRACK LENGTH WAS: " << trueLength << std::endl;
  std::cout << "*************************************************" << std::endl;

  path.SetTrueTrackLength(trueLength);
}


void trex::TTRExPIDAlgorithm::TrackCleanliness(trex::TTRExPath& path){
  
  std::vector<int> TrackNumbers;
  std::vector<int> TrackFractions;
  std::vector<trex::TTrueTrack*> TrueTracks;
  int TotalNumberOfPathHits=0;

  std::vector<trex::TTRExHVCluster*>& Clusters = path.GetClusters();
  for(auto clusterit=Clusters.begin(); clusterit != Clusters.end(); clusterit++){
    
    trex::TTRExHVCluster& Cluster=**clusterit;
    std::vector<trex::TTPCHitPad*>& Hits = Cluster.GetClusterHits();

    for(auto hitit=Hits.begin(); hitit != Hits.end(); hitit++){
      
      trex::TTPCHitPad& Hit = **hitit;
     
      TotalNumberOfPathHits+=1;

      trex::TTrueTrack* trueTrack = Hit.GetTrueTrack();
      if(!trueTrack) continue;
      int TrackNumber = trueTrack->GetTrackNumber();
      bool found=false;

      if(TrackNumbers.size()==0){
	TrackNumbers.push_back(TrackNumber);
	TrackFractions.push_back(1);
	TrueTracks.push_back(trueTrack);
      }
      else{
	for(int i=0; i<TrackNumbers.size(); ++i){
	  if (TrackNumber==TrackNumbers[i]){
	    TrackFractions[i]+=1;
	    found=true;
	  }
	}
	if(!found){
	  
	  std::cout << "TrackNumber was not found, create new entry" << std::endl; 
	  TrackNumbers.push_back(TrackNumber);
	  TrackFractions.push_back(1);
	  TrueTracks.push_back(Hit.GetTrueTrack());
	}
      }
    }
  }
  if(!TrackFractions.size()){
    path.SetNumberOfTrueHitsFound(0);
    path.SetTrackCompleteness(0);
    path.SetTrackCleanliness(0);
    return;
  }

  int largestFraction=0;
  int mainTrackIndex=0;
 
  for(int j=0; j<TrackFractions.size(); ++j){
    
    if(largestFraction < TrackFractions[j]){      
      std::cout << "Track Fraction " << j << " is " << TrackFractions[j] << std::endl;
      largestFraction=TrackFractions[j];
      mainTrackIndex = j;
    }
  }

  path.FillFromTruthTrack(TrueTracks[mainTrackIndex]);
  path.SetNumberOfTrueHitsFound(largestFraction);
  
  int TrueNumberOfHits = TrueTracks[mainTrackIndex]->GetNumberOfHits();
 
  double completeness = (double)largestFraction/(double)TrueNumberOfHits;

  std::cout << "The completeness was: " << completeness << std::endl;

  if(TotalNumberOfPathHits < largestFraction){
    std::cout << "THIS IS NOT RIGHT!" << std::endl;
    std::cout << "Hits in Path is: " << TotalNumberOfPathHits << " and largest Fraction is: " << largestFraction << std::endl; 
  }
  
  double cleanliness = (double)largestFraction/(double)TotalNumberOfPathHits; 
  
  path.SetTrackCompleteness(completeness);
  path.SetTrackCleanliness(cleanliness);

}


//Set whether it is pro or pi according to given cut (value defined in PIDAlgo.hxx
void trex::TTRExPIDAlgorithm::PID(trex::TTRExPath& path){
  
  double fPIDcut = 3.5e-6;
  double dEdx = path.GetdEdx();
  
  if(dEdx > fPIDcut){
    path.SetPID(-1);  
    std::cout << "This is a proton" << std::endl;
  }else{path.SetPID(1);
    std::cout << "This is a pion" << std::endl;}
  
}

