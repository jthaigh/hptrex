
// eddy
#include "TEventDisplay.hxx"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TDirectory.h"
#include<map>

trex::TEventDisplay::TEventDisplay(TFile* plotFile,unsigned int nEvt) {

  fPlotFile=plotFile;
  iEvt=nEvt;
}

void trex::TEventDisplay::Process(std::vector<trex::TTPCHitPad*>& hits, std::vector<TTrueHit*>& trueHits, trex::TTRExEvent* event, trex::TTPCLayout& layout){
  
  std::vector<trex::TTRExPattern>& patterns=event->GetPatterns();

  std::vector<TGraph*> xyGraphs;
  std::vector<TGraph*> xzGraphs;

  // get patterns

  gStyle->SetOptStat(0);
  int iColor=0;
  int colors[11]={kBlue, kYellow, kGreen, kMagenta, kCyan, kOrange, kPink, kAzure, kSpring, kViolet};
    
  // set up container for hitpad level unused  
  std::vector<trex::TTPCHitPad*> usedTREx;
  
  for(auto patIt = patterns.begin(); patIt != patterns.end(); ++patIt){
    trex::TTRExPattern& pat = *patIt;
    std::vector<trex::TTRExPath>& subPaths= pat.GetPaths();
    std::vector<trex::TTRExJunction>& subJuncts=pat.GetJunctions();

    //PD NEED TO PUT BETTER FILLING METHOD HERE
    //extract the path hits from the sub algorithm and build HVCluster objects to fill the paths with

    for(auto iPath=subPaths.begin(); iPath!=subPaths.end(); ++iPath) {
      
      std::vector<trex::TTRExHVCluster*> clusters = iPath->GetClusters();
      for(auto iCluster=clusters.begin(); iCluster!=clusters.end(); ++iCluster){

	vector<trex::TTPCHitPad*> cHits = (*iCluster)->GetClusterHits();

	for(auto iHit=cHits.begin(); iHit!=cHits.end(); ++iHit){

	  if(std::find(usedTREx.begin(), usedTREx.end(), *iHit)==usedTREx.end()){
	    usedTREx.push_back(*iHit);
	  }
	}	
      }
    }

    //PD NEED TO PUT BETTER FILLING METHOD HERE
    //Fill the junctsContainer with junctions fromt he subevent
    
    for(auto iJunct=subJuncts.begin(); iJunct!=subJuncts.end(); ++iJunct) {      
            
      std::vector<trex::TTPCHitPad*> hits = iJunct->GetHits();

      for(auto iHit=hits.begin(); iHit!=hits.end(); ++iHit) {
	
	//fill usedTREx
	if(std::find(usedTREx.begin(),usedTREx.end(),*iHit)==usedTREx.end()){
	  usedTREx.push_back(*iHit);
	}
      }
    }
    

    for(auto iPath=subPaths.begin();iPath!=subPaths.end();++iPath){
      
      std::vector<trex::TTPCHitPad*> usedThisObject;
      int color_index = iColor%10;
      int color_increment = iColor%4;
      xyGraphs.push_back(new TGraph(1));
      xzGraphs.push_back(new TGraph(1));
      xyGraphs.back()->SetMarkerColor(colors[color_index]+color_increment);
      xyGraphs.back()->SetMarkerStyle(20);
      xyGraphs.back()->SetMarkerSize(0.5);
      xzGraphs.back()->SetMarkerColor(colors[color_index]+color_increment);
      xzGraphs.back()->SetMarkerStyle(20);
      xzGraphs.back()->SetMarkerSize(0.5);
      iColor++;
      unsigned int iPt=0;

      for(auto iClu=iPath->GetClusters().begin();iClu!=iPath->GetClusters().end();++iClu){
	std::vector<trex::TTPCHitPad*> hits=(*iClu)->GetClusterHits();
	for(auto iHit=hits.begin();iHit!=hits.end();++iHit){
	  TVector3 pos=(*iHit)->GetPosition();
	  xyGraphs.back()->SetPoint(iPt,pos.Y(),pos.Z());
	  xzGraphs.back()->SetPoint(iPt++,pos.X(),pos.Z());
	  
	}
      }
    }
    
    for(auto iJunct=subJuncts.begin();iJunct!=subJuncts.end();++iJunct){
      
      //int color_index = iColor%10;
      //int color_increment = iColor%4;
      xyGraphs.push_back(new TGraph(1));
      xzGraphs.push_back(new TGraph(1));
      xyGraphs.back()->SetMarkerColor(kRed);
      xyGraphs.back()->SetMarkerStyle(21);
      xyGraphs.back()->SetMarkerSize(0.5);
      xzGraphs.back()->SetMarkerColor(kRed);
      xzGraphs.back()->SetMarkerStyle(21);
      xzGraphs.back()->SetMarkerSize(0.5);
      //iColor++;
      unsigned int iPt=0;
      std::vector<trex::TTPCHitPad*> hits=iJunct->GetHits();
      for(auto iHit=hits.begin();iHit!=hits.end();++iHit){
	TVector3 pos=(*iHit)->GetPosition();
	xyGraphs.back()->SetPoint(iPt,pos.Y(),pos.Z());
	xzGraphs.back()->SetPoint(iPt++,pos.X(),pos.Z());
      }
    }
  }

  std::vector<int> strangePDG;
  
  std::map<int,int> truthColors;
    
  truthColors[13]=kSpring-9;
  truthColors[-13]=kOrange-3;
  truthColors[11]=kPink-4;
  truthColors[-11]=kPink-7;
  truthColors[22]=kOrange-2;
  truthColors[211]=kCyan-9;
  truthColors[-211]=kAzure-3;
  truthColors[2112]=kYellow-4;
  truthColors[2212]=kRed-4;
  
  int trueTrackCount=0;
  
  std::map<int,TGraph*> xyHitGraphs;
  std::map<int,TGraph*> xzHitGraphs;

  for (auto iTrueHits=trueHits.begin(); iTrueHits!=trueHits.end();++iTrueHits){
  
  int trackId = (*iTrueHits)->TrueTrackID;
  
  if(!xyHitGraphs.count(trackId)){
    
    xyGraphs.push_back(new TGraph(0));
    xzGraphs.push_back(new TGraph(0));
    xyHitGraphs[trackId]=xyGraphs.back();
    xzHitGraphs[trackId]=xzGraphs.back();

    int pdg = (*iTrueHits)->pdg;
    int color=0;

    if(truthColors.find(pdg) == truthColors.end()){
      strangePDG.push_back(pdg);
    }
    else{color = truthColors[pdg];}

    xyGraphs.back()->SetMarkerColor(color);
    xyGraphs.back()->SetMarkerStyle(31);
    xyGraphs.back()->SetMarkerSize(0.1);
    xzGraphs.back()->SetMarkerColor(color);
    xzGraphs.back()->SetMarkerStyle(31);
    xzGraphs.back()->SetMarkerSize(0.1);

    trueTrackCount++;
  }
    
  TLorentzVector pos = (*iTrueHits)->TruePos4;
  
  //xyHitGraphs[trackId]->SetPoint(xyHitGraphs[trackId]->GetN(),0.1*pos.X(),0.1*pos.Y());
  //xzHitGraphs[trackId]->SetPoint(xzHitGraphs[trackId]->GetN(),0.1*pos.X(),0.1*pos.Z());
  //0.1* factor turns positions from mm to cm
  
  xyHitGraphs[trackId]->SetPoint(xyHitGraphs[trackId]->GetN(),pos.Y(),pos.Z());
  xzHitGraphs[trackId]->SetPoint(xzHitGraphs[trackId]->GetN(),pos.X(),pos.Z());

  }
  
  if(strangePDG.size()!=0){
    std::cout << "STRANGE PDGS FOUND!!!" << std::endl;
    for(auto istrange=strangePDG.begin(); istrange!=strangePDG.end();++istrange){
      std::cout << "STRANGE: " << *istrange << std::endl;}
  }

 

  if(hits.size()){
    fPlotFile->cd();
    
    std::vector<trex::TTPCHitPad*> unusedHits;  
    for(auto iHit=hits.begin();iHit!=hits.end();++iHit){
      if(std::find(usedTREx.begin(),usedTREx.end(),*iHit)==usedTREx.end()){
	unusedHits.push_back(*iHit);
      }
    }

    xyGraphs.push_back(new TGraph(1));
    xzGraphs.push_back(new TGraph(1));
    xyGraphs.back()->SetMarkerColor(1);
    xyGraphs.back()->SetMarkerStyle(20);
    xyGraphs.back()->SetMarkerSize(0.2);
    xzGraphs.back()->SetMarkerColor(1);
    xzGraphs.back()->SetMarkerStyle(20);
    xzGraphs.back()->SetMarkerSize(0.2);
    
    unsigned int iPt=0;
    
    for(auto iHit=unusedHits.begin();iHit!=unusedHits.end();++iHit){
      TVector3 pos=(*iHit)->GetPosition();
      xyGraphs.back()->SetPoint(iPt,pos.Y(),pos.Z());
      xzGraphs.back()->SetPoint(iPt++,pos.X(),pos.Z());
    }

    char buf[20];
    sprintf(buf,"evt_%d_xy",iEvt);
    TCanvas cxy(buf,buf);
        
    TH2F dummyxy("XY-view","XY-view",
		 1000,layout.GetMinPos().Y()-10.,layout.GetMaxPos().Y()+10.,
		 1000,layout.GetMinPos().Z()-10.,layout.GetMaxPos().Z()+10.);
    
    dummyxy.Draw();

    for(auto iGr=xyGraphs.begin();iGr!=xyGraphs.end();++iGr){
      cxy.cd();
      (*iGr)->Draw("Psame");
      double x,y;
      (*iGr)->GetPoint(0,x,y);
    }
    
    cxy.Write();

    sprintf(buf,"evt_%d_xz",iEvt);
    TCanvas cxz(buf,buf);
    
    TH2F dummyxz("XZ-view","XZ-view",
		 1000,layout.GetMinPos().X()-10.,layout.GetMaxPos().X()+10.,
		 1000,layout.GetMinPos().Z()-10.,layout.GetMaxPos().Z()+10.);
    
    
    dummyxz.Draw();

    for(auto iGr=xzGraphs.begin();iGr!=xzGraphs.end();++iGr){
      cxz.cd();
      (*iGr)->Draw("Psame");
      double x,y;
      (*iGr)->GetPoint(0,x,y);
    }

    cxz.Write();

    for(auto iGr =xyGraphs.begin();iGr!=xyGraphs.end();++iGr) delete *iGr;
    for(auto iGr=xzGraphs.begin();iGr!=xzGraphs.end();++iGr) delete *iGr;
  }
  
}
