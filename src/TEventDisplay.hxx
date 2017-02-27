#ifndef TEVENTDISPLAY_HXX
#define TEVENTDISPLAY_HXX

// c++
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>

// ROOT
#include <TFile.h>
#include <TStyle.h>

// TREx
#include "TTrueHit.hxx"
#include "TTRExPattern.hxx"
#include "TTRExPath.hxx"
#include "TTRExJunction.hxx"
#include "TTPCLayout.hxx"

namespace trex{
  /// Pattern recognition for TPC tracks based on path finding.
  class TEventDisplay {
    public:
      /// Default constructor
    TEventDisplay(TFile* plotFile,unsigned int nEvt);
    /// Clean up values from previous processing
    
    void Process(std::vector<trex::TTPCHitPad*>& hits, std::vector<TTrueHit*>& trueHits, trex::TTRExEvent* event, trex::TTPCLayout& layout);

    private:

    TFile* fPlotFile;
    unsigned int iEvt;
    
  };
}

#endif
