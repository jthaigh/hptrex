//------------------------------------------------------------------

// ALL NEW DATA CLASS TYPES FOR OUTPUT PATTERN

//------------------------------------------------------------------





#ifdef __CINT__



//Output Classes

                                                           
#pragma link C++ class trex::TTRExEvent+;                                  
#pragma link C++ class trex::TTRExPattern+;                           
#pragma link C++ class std::vector<trex::TTRExPattern>+; 
                  
#pragma link C++ class trex::TTPCHitPad+;

#pragma link C++ class trex::TTPCHitPad*+;
 
#pragma link C++ class std::vector<trex::TTPCHitPad>+; 

#pragma link C++ class trex::TTRExHVCluster+;
#pragma link C++ class std::vector<trex::TTRExHVCluster>+;
#pragma link C++ class std::vector<std::vector<trex::TTRExHVCluster> >+;

#pragma link C++ class trex::TTRExPath+;
#pragma link C++ class std::vector<trex::TTRExPath>+;

                   
//#pragma link C++ class std::vector<std::vector<trex::TTPCHitPad*> >+;


#endif
