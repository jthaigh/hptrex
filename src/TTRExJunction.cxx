#include "TTRExJunction.hxx"


void trex::TTRExJunction::SetConnectedPaths(std::vector<trex::TTRExPath*> &paths) {
  
  fConnectedPaths.clear();
  fConnectedPathsId.clear();
  
  fConnectedPaths = paths;
  for(int i=0; i<paths.size(); ++i){
    fConnectedPathsId.push_back(paths[i]->GetId());  
  }
}


void trex::TTRExJunction::AddConnectedPath(trex::TTRExPath* path){                                                                                    
  if(find(fConnectedPaths.begin(), fConnectedPaths.end(), path) != fConnectedPaths.end()){                                              
    std::cout << "This Path is already connected" << std::endl;                                                                         
    return;                                                                                                                             
  }                                                                                                                                     
  fConnectedPaths.push_back(path);                                                                                                      
  fConnectedPathsId.push_back(path->GetId());                                                                                           
} 


bool trex::TTRExJunction::IsPathConnected(unsigned int WantedPathId){
  bool isConnected = false;
  for(int i=0; i<fConnectedPathsId.size(); ++i){
    if(fConnectedPathsId[i] == WantedPathId){
      isConnected = true;
      break;
    }
  }
  return isConnected;
}

void trex::TTRExJunction::SetHits(const std::vector<TTPCHitPad*>& theHits){
  fHits=theHits;
  TVector3 position;
  double totalCharge=0;
  for(auto iHit=fHits.begin();iHit!=fHits.end();++iHit){
    position+=(*iHit)->GetPosition()*(*iHit)->GetCharge();
    totalCharge+=(*iHit)->GetCharge();
  }
  fPosition=(1./totalCharge)*position;
}
