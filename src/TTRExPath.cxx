#include "TTRExPath.hxx"



void trex::TTRExPath::SetConnectedJunctions(std::vector<trex::TTRExJunction*> &juncts){
  fConnectedJunctions.clear();
  fConnectedJunctionsId.clear();
  
  fConnectedJunctions = juncts;
  for(int i=0; i<juncts.size(); ++i){
    fConnectedJunctionsId.push_back(juncts[i]->GetId());
  }
}


void trex::TTRExPath::AddConnectedJunction(trex::TTRExJunction* junct){

  if(find(fConnectedJunctions.begin(),fConnectedJunctions.end(), junct) != fConnectedJunctions.end()){
    std::cout << "This Junction is already connected" << std::endl;
    return;
  }
  fConnectedJunctions.push_back(junct);
  fConnectedJunctionsId.push_back(junct->GetId());                                                                                    
}
