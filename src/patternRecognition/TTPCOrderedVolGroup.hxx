#ifndef TTPCORDEREDVOLGROUP_HXX
#define TTPCORDEREDVOLGROUP_HXX

// c++
#include <set>
#include <map>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <utility>
#include <algorithm>

// ROOT
#include <TObject.h>
#include <TMath.h>

// nd280
#include <THandle.hxx>

// TREx
#include <TTPCHVCluster.hxx>
#include <TTPCDebug.hxx>

// eddy
#include "TTPCLayout.hxx"
#include "TTPCPathVolume.hxx"
#include "TTPCUnitVolume.hxx"
#include "TTPCVolGroup.hxx"

namespace ND{
  class TTPCOrderedVolGroup;
}

/// class for holding an ordered chain of cells between two groups, and also holding pointers to the groups themselves
class ND::TTPCOrderedVolGroup : public TObject {
  public:
    /// Default constructor
    TTPCOrderedVolGroup();
    /// Constructor
    TTPCOrderedVolGroup(ND::TTPCLayout* fLayout);
    /// Destructor
    virtual ~TTPCOrderedVolGroup();

    /// Add a pointer to the hits corresponding to hits near the start of this group (hits near hit at size()-1 index)
    void AddFrontHits(ND::THandle<ND::TTPCVolGroup> frontHits);
    /// Add a pointer to the hits corresponding to hits near the end of this group (hits near hit at 0 index)
    void AddBackHits(ND::THandle<ND::TTPCVolGroup> backHits);
    /// Add a pointer to extended group of hits associated with this path
    void AddExtendedHits(ND::THandle<ND::TTPCVolGroup> extendedHits);

    /// Check front hits exist and are not empty
    bool HasFrontHits(){ return fAddedFrontHits ? !fFrontHits->empty() : false; }
    /// Check back hits exist and are not empty
    bool HasBackHits(){ return fAddedBackHits ? !fBackHits->empty() : false; }
    /// Get average position of front hits 
    TVector3 GetFrontAveragePosition(){ return fFrontHits->GetAveragePosition(); }
    /// Get average position of back hits
    TVector3 GetBackAveragePosition(){ return fBackHits->GetAveragePosition(); }
    /// Get average cell in front hits
    ND::TTPCUnitVolume* GetFrontAverageVol(){ return fFrontHits->GetAverageVol(); }
    /// Get average cell in back hits
    ND::TTPCUnitVolume* GetBackAverageVol(){ return fBackHits->GetAverageVol(); }
    /// Return bare pointer to hits in front
    ND::THitSelection* GetFrontHitSelection(){ return fFrontHits->GetHitSelection(); }
    /// Return bare pointer to hits in back
    ND::THitSelection* GetBackHitSelection(){ return fBackHits->GetHitSelection(); }
    /// Get handle for front hits
    ND::THandle<ND::TTPCVolGroup> GetFrontHits(){ return fFrontHits; }
    /// Get handle for back hits
    ND::THandle<ND::TTPCVolGroup> GetBackHits(){ return fBackHits; }
    /// Get ID for the group of hits at the front of this group
    unsigned int GetFrontID(){ return fFrontHits->GetID(); }
    /// Get ID for the group of hits at the back of this group
    unsigned int GetBackID(){ return fBackHits->GetID(); }

    /// Get extended group of hits associated with this path
    //MDH
    //Changed this to copy fExtendedHits into input object.
    //Slightly ugly but I think this may be the right way
    //since we want to manage persistence separately.
    void GetExtendedHits(ND::TTPCVolGroup& ret){  ret=fExtendedHits; }
  
    /// Iterator to start of front hits
    std::map<long, ND::TTPCUnitVolume*>::iterator frontHitsBegin(){ return fFrontHits->begin(); }
    /// Iterator to end of front hits
    std::map<long, ND::TTPCUnitVolume*>::iterator frontHitsEnd(){ return fFrontHits->end(); }
    /// Iterator to start of back hits
    std::map<long, ND::TTPCUnitVolume*>::iterator backHitsBegin(){ return fBackHits->begin(); }
    /// Iterator to end of back hits
    std::map<long, ND::TTPCUnitVolume*>::iterator backHitsEnd(){ return fBackHits->end(); }
    /// Iterator to start of extended hits
    std::map<long, ND::TTPCUnitVolume*>::iterator extendedHitsBegin(){ return fExtendedHits->begin(); }
    /// Iterator to end of extended hits
    std::map<long, ND::TTPCUnitVolume*>::iterator extendedHitsEnd(){ return fExtendedHits->end(); }

    /// Push a new path volume back into the vector of hits in this group and return it
    ND::TTPCPathVolume* AddCell(ND::TTPCUnitVolume* cell, bool isXCluster=false);
    /// Associate hits from input group with hits from the main path
    void DoClustering(bool partial=false);
    /// General clustering for most paths 
    void DoStandardClustering(bool partial=false);
    /// General clustering for most paths, ensuring that all hits that can be merged in are merged in
    void DoGreedyClustering(bool partial=false);
    /// Specialised clustering for paths long in x
    void DoXClustering();
    /// Remove sequences of hits changing only in x to simplify HV clustering procedure
    void StripX();
    /// Use dichotomy technique to determine local angle at each point in a cluster
    void FindClusterAnglesByDichotomy();
    /// Turn local cluster angles into classification as horizontal or vertical clusters
    void FindHVClustersFromAngles();
    /// Extrapolate HV clusters to include skipped or otherwise missed hits
    void ExtrapolateHVClusters();
    /// Interpolate HV clusters to include skipped or otherwise missed hits
    void InterpolateHVClusters();

    /// Fill clusters depending on whether they're classed as horizontal or vertical
    void FillHVClusters();
    /// Merge hits from main path that are in same direction and overlap
    void MergeHVClusters();

    /// Greedily fill clusters depending on whether they're classed as horizontal or vertical
    void GreedyFillHVClusters();
    /// Merge directly adjacent clusters for greedy filling
    void GreedyMergeHVClusters(std::map<int, std::vector<ND::TTPCPathVolume*> >& clusters);
    /// Find path volume with most average position in a set
    ND::TTPCPathVolume* GetAverageVol(std::vector<ND::TTPCPathVolume*> clusters);
    /// Try to add a volume to a cluster
    bool GreedyAddVol(ND::TTPCUnitVolume* vol, std::map<int, std::vector<ND::TTPCPathVolume*> >& clusters, bool isVertical);

    /// Create and populate new x clusters
    void ExpandXClusters();
    /// Set x cluster angles to default
    void SetXClusterHVAngles();

    /// Clean up empty or null cells
    void Clean();
    /// Reverse the order of elements in the path
    void Flip();
    /// Order in forwards track direction
    void OrderForwardsDirection();
    /// Order for negative curvature
    void OrderNegativeCurvature();
    /// Order from junction position
    void OrderFromJunction(ND::THandle<ND::TTPCVolGroup> junction);
    /// Order from position
    void OrderFromPosition(TVector3 pos);

    /// Get whether the path meets criteria for delta ray tagging
    bool GetDeltaCriteriaMet();
    /// Get list of all hits associated with elements of this group
    ND::THandle<ND::THitSelection> GetClusters();

    /// Set front hits to be a vertex (or not)
    void SetFrontIsVertex(bool isVertex=false){ fFrontIsVertex = isVertex; }
    /// Set back hits to be a vertex (or not)
    void SetBackIsVertex(bool isVertex=false){ fBackIsVertex = isVertex; }
    /// Get status of front hits as vertex
    bool GetFrontIsVertex() { return fFrontIsVertex; }
    /// Get status of back hits as vertex
    bool GetBackIsVertex() { return fBackIsVertex; }
    /// Get whether the group has extended hits associated with it or not
    bool GetHasExtendedHits(){ return fHasExtendedHits; }

    /// Get minimum x in path hits
    int GetXMin(){Close(); return fXMin;}
    /// Get maximum x in path hits
    int GetXMax(){Close(); return fXMax;}
    /// Get minimum y in path hits
    int GetYMin(){Close(); return fYMin;}
    /// Get maximum y in path hits
    int GetYMax(){Close(); return fYMax;}
    /// Get minimum z in path hits
    int GetZMin(){Close(); return fZMin;}
    /// Get maximum z in path hits
    int GetZMax(){Close(); return fZMax;}
    /// Get x extent of path hits
    int GetXSize(){Close(); return (1+fXMax-fXMin);}
    /// Get y extent of path hits
    int GetYSize(){Close(); return (1+fYMax-fYMin);}
    /// Get z extent of path hits
    int GetZSize(){Close(); return (1+fZMax-fZMin);}
    /// Get whether this cluster counts as an x path
    bool GetIsXPath(){ return fIsXPath; }

    /// Set whether this cluster counts as an x path
    void SetIsXPath(bool isXPath){ fIsXPath = isXPath; }

    /// Erase element from hits in this group
    void erase(std::vector<ND::TTPCPathVolume*>::iterator hit){ fHits.erase(hit); fClosed = false; }
    /// Erase elements from hits in this group
    void erase(std::vector<ND::TTPCPathVolume*>::iterator begin, std::vector<ND::TTPCPathVolume*>::iterator end){ fHits.erase(begin, end); fClosed = false; }
    /// Returns whether hit vector is empty 
    bool empty(){ return fHits.empty(); }
    /// Returns size of hit vector
    int size(){ return int(fHits.size()); }
    /// Returns element at a specified index
    ND::TTPCPathVolume* at(int i){ return fHits.at(i); }
    /// Return element at i
    ND::TTPCPathVolume*& operator[](int i){ return fHits[i]; }
    /// Iterator to begining of hit map
    std::vector<ND::TTPCPathVolume*>::iterator begin(){ return fHits.begin(); }
    /// Iterator to end of hit map
    std::vector<ND::TTPCPathVolume*>::iterator end(){ return fHits.end(); }
    /// Reverse iterator to end of hit map
    std::vector<ND::TTPCPathVolume*>::reverse_iterator rbegin(){ return fHits.rbegin(); }
    /// Reverse iterator to beginning of hit map
    std::vector<ND::TTPCPathVolume*>::reverse_iterator rend(){ return fHits.rend(); }

    /// Get string representing cluster orientations
    std::string GetOrientations();
    /// Print cluster orientations
    void PrintOrientations();
    /// Print cluster positions
    void PrintPositions(bool size=false, bool clusterHits=false);
    /// Debugging tool for checking hits for memory issues
    void CheckClusters();

  private:
    /// TPC layout associated with this group
    ND::TTPCLayout* fLayout;

    /// List of hits in this group
    std::vector<ND::TTPCPathVolume*> fHits;

    /// Hits near the front of this group (near hit at size()-1 index)
    ND::THandle<ND::TTPCVolGroup> fFrontHits;
    /// Hits near the end of this group (near hit at 0 index)
    ND::THandle<ND::TTPCVolGroup> fBackHits;

    /// Hits associated with the whole path
    ND::THandle<ND::TTPCVolGroup> fExtendedHits;

    /// Whether front hits have been added
    bool fAddedFrontHits;
    /// Whether back hits have been added
    bool fAddedBackHits;

    /// Whether front hits represent a vertex
    bool fFrontIsVertex;
    /// Whether back hits represent a vertex
    bool fBackIsVertex;

    /// Whether the group has extended hits associated with it or not
    bool fHasExtendedHits;

    /// Minimum x in path hits
    int fXMin;
    /// Maximum x in path hits
    int fXMax;
    /// Minimum y in path hits
    int fYMin;
    /// Maximum y in path hits
    int fYMax;
    /// Minimum z in path hits
    int fZMin;
    /// Maximum z in path hits
    int fZMax;

    /// Whether this cluster counts as an x cluster
    bool fIsXPath;
    /// Whether internal variables have been calculated (open when fHits changes)
    bool fClosed;

    /// Set up initial variables
    void SetUp();
    /// Recursive function to use for dichotomy technique of working out angle
    void RecursiveDichotomy(std::vector<float>& prevAngles, int firstID, int lastID, float prevAngDiff);
    /// Extrapolate HV clusters in a specific dirction
    std::vector<ND::TTPCPathVolume*> GetExtrapolatedClusters(std::vector<ND::TTPCPathVolume*>::iterator pathVolIt, int dir);

    /// Calculate internal variables
    void Close();

    /// Export class for ROOT
    ClassDef(TTPCOrderedVolGroup, 1);
};

#endif
