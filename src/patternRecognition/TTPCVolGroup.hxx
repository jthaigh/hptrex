#ifndef TTPCVOLGROUP_HXX
#define TTPCVOLGROUP_HXX

// c++
#include <map>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

// ROOT
#include <TObject.h>
#include <TVector3.h>

#include <cstdlib>

// eddy
#include "TTPCLayout.hxx"
#include "TTPCUnitVolume.hxx"

namespace trex{
  class TTPCVolGroup;
}

/// class for holding and operating on groups of the cells used in feature finding and path finding
class trex::TTPCVolGroup{
  public:
    /// Constructor
    TTPCVolGroup(trex::TTPCLayout* layout, unsigned int id=0);
    /// Destructor
    virtual ~TTPCVolGroup();

    /// Add whole map of hits to group from scratch
    void AddHitMap(std::map<long, trex::TTPCUnitVolume*> hitMap);
    /// Copy hits in an existing map into this group's map
    void AddHits(std::map<long, trex::TTPCUnitVolume*> hitMap);
    /// Copy hits between two iterators into this group's map 
    void AddHits(std::map<long, trex::TTPCUnitVolume*>::iterator hitMapBegin, std::map<long, trex::TTPCUnitVolume*>::iterator hitMapEnd);
    /// Copy hits in from a vector into this group's map
    void AddHits(std::vector<trex::TTPCUnitVolume*> hitList);
    /// Copy hits between two vector iterators into this group's map 
    void AddHits(std::vector<trex::TTPCUnitVolume*>::iterator hitListBegin, std::vector<trex::TTPCUnitVolume*>::iterator hitListEnd);
    /// Copy hits in another group's map into this group's map
    void AddHits(trex::TTPCVolGroup& hits);
    /// Add individual hit, returning false if it already existed in this group's map
    bool AddHit(trex::TTPCUnitVolume* hit, bool safe=true);

    /// Remove a hit by its id from this group's map, returning false if it wasn't found and true if it was sucessfully removed
    bool RemoveHit(long id);
    /// Mark a hit for deletion by setting its pointer to zero
    void MarkHit(long id);
    /// Clear hits marked for deletion by setting their pointers to zero
    void ClearMarked();

    /// Merge one group into this group, summing charges where cells overlap
    void MergeHits(trex::TTPCVolGroup& hits);

  //MDH
  //Not used
  /*
    /// Add 'fake' cell not corresponding to a detected hit
    void AddPseudoHit(trex::TTPCUnitVolume* hit);
    /// Add 'fake' cell not corresponding to a detected hit from its x, y and z id and charge
    void AddNewPseudoHit(int x, int y, int z, float q);
    /// Add gaussian distribution of 'fake' cells in two dimensions (x view, y view or z view with axis 1, 2 or 3 respectively), around specific x, y and z id with specific x, y and z spread
    void AddNewPseudoGauss(int x, int y, int z, float q, int axis, float sigmaX, float sigmaY, float sigmaZ);
  */
    /// Find nearest hit in group to a given x, y and z
    long GetNearestHit(int x, int y, int z, int maxDist=1000);
    /// Find position of nearest hit in group to a given x, y and z
    TVector3 GetNearestHitPos(int x, int y, int z, int maxDist);

    /// Get list of all hits associated with elements of this group
  std::vector< trex::TTPCHitPad* > GetHits();

    /// Get ID corresponding to the track end or vertex this junction is supposed to represent
    unsigned int GetID(){ return fID; }

    /// Return whether this group contains an element corresponding to the supplied unique id
    bool Contains(long id);
    /// Return whether this group contains an element corresponding to the supplied volume element
    bool Contains(TTPCUnitVolume* vol);
    /// Get wheter a hit is at a lean in this group
    bool GetLeanValid(trex::TTPCUnitVolume* vol);

    /// Return an element from hit map 
    trex::TTPCUnitVolume* GetHit(long id, bool safe=true);
    /// Return the hit map
    std::map<long, trex::TTPCUnitVolume*>& GetHitMap(){ return fHitMap; }

    /// Return average position of all cells in the group 
    TVector3 GetAveragePosition(){ Close(); return fAveragePosition; }
    /// Get vector representing average cell x, y, z
    TVector3 GetAveragePosXYZ(){ Close(); return TVector3( (double)fAveragePad.x, (double)fAveragePad.y, (double)fAveragePad.z ); }
    /// Get the pad at the average position of all cells in the group
    trex::TTPCCellInfo3D GetAveragePad(){ Close(); return fAveragePad; }
    /// Get the unit volume at the average position of all cells in the group
    trex::TTPCUnitVolume* GetAverageVol(){ Close(); return fAverageUnitVolume; }

  //MDH
  //Redundant under new way of storing groups of hits
    /// Get bare pointer to the hit selection containing all cells in the group
  //    trex::THitSelection* GetHitSelection();
    /// Get the average cell x
    int GetX(){ Close(); return fAveragePad.x; }
    /// Get the average cell y
    int GetY(){ Close(); return fAveragePad.y; }
    /// Get the average cell z
    int GetZ(){ Close(); return fAveragePad.z; }
    /// Get standard deviation of all pads in this group along x axis
    float GetSigmaPadX(){ Close(); return fSigmaPadX; }
    /// Get standard deviation of all pads in this group along y axis
    float GetSigmaPadY(){ Close(); return fSigmaPadY; }
    /// Get standard deviation of all pads in this group along z axis
    float GetSigmaPadZ(){ Close(); return fSigmaPadZ; }
    /// Get the unique id of the pad at the average position in this group
    long GetAveragePadID(){ Close(); return fAveragePadID; }
    /// Get the total charge contained in this group
    float GetCharge(){ Close(); return fCharge; }
    /// Get average charge contained in this group's cells
    float GetAverageCharge(){ return fAverageCharge; }

    /// Get filter on x position to use for calculating average cell
    int GetXLean(){ Close(); return fXLean; }
    /// Get filter on y position to use for calculating average cell
    int GetYLean(){ Close(); return fYLean; }
    /// Get filter on z position to use for calculating average cell
    int GetZLean(){ Close(); return fZLean; }

    /// Get min pad x in this group
    int GetXMin(){ Close(); return fXMin; }
    /// Get max pad x in this group
    int GetXMax(){ Close(); return fXMax; }
    /// Get difference between min and max pad x in this group
    int GetXSize(){ Close(); return fXSize; }
    /// Get min pad y in this group
    int GetYMin(){ Close(); return fYMin; }
    /// Get max pad y in this group
    int GetYMax(){ Close(); return fYMax; }
    /// Get difference between min and max pad y in this group
    int GetYSize(){ Close(); return fYSize; }
    /// Get min pad z in this group
    int GetZMin(){ Close(); return fZMin; }
    /// Get max pad z in this group
    int GetZMax(){ Close(); return fZMax; }
    /// Get difference between min and max pad z in this group
    int GetZSize(){ Close(); return fZSize; }

    /// Returns size of hit map
    int size(){ return int(fHitMap.size()); }
    /// Check if this group's map has size zero
    bool empty(){ return fHitMap.empty(); }
    /// Clear this group's map
    void clear(){ fIsClosed = false; fHitMap.clear(); }
    /// Iterator to begining of hit map
    std::map<long, trex::TTPCUnitVolume*>::iterator begin(){ return fHitMap.begin(); }
    /// Iterator to end of hit map
    std::map<long, trex::TTPCUnitVolume*>::iterator end(){ return fHitMap.end(); }
    /// Find an element
    std::map<long, trex::TTPCUnitVolume*>::iterator find(long id){ return fHitMap.find(id); }
    /// Return an element
    trex::TTPCUnitVolume* GetEl(long id){ std::map<long, trex::TTPCUnitVolume*>::iterator vol = find(id); return ((vol == end()) ? 0 : vol->second); }

    /// Set status for filtering out x hits when getting positions
    void SetXLean(int xLean=0, bool force=false){ SetLean(fXLean, xLean, force); }
    /// Set status for filtering out y hits when getting positions
    void SetYLean(int yLean=0, bool force=false){ SetLean(fYLean, yLean, force); }
    /// Set status for filtering out z hits when getting positions
    void SetZLean(int zLean=0, bool force=false){ SetLean(fZLean, zLean, force); }
    /// General function for setting a lean value
    void SetLean(int& leanVar, int leanVal=0, bool force=false);

    /// Check for any empty or broken hits
    void CheckHits();

    /// Static function for returning a free ID
    static unsigned int GetFreeID(){ return ++sMaxID; }
    /// Static function for reseting free ID
    static void ResetFreeID(){ sMaxID = 0; }

  private:
    /// Layout in which the cells in this group live
    trex::TTPCLayout* fLayout;

    /// Filter on x position to use for calculating average cell
    int fXLean;
    /// Filter on y position to use for calculating average cell
    int fYLean;
    /// Filter on z position to use for calculating average cell
    int fZLean;

    /// Minimum x value in group
    int fXMin;
    /// Maximum x value in group
    int fXMax;
    /// Difference between min and max x values in group
    int fXSize;
    /// Minimum y value in group
    int fYMin;
    /// Maximum y value in group
    int fYMax;
    /// Difference between min and max y values in group
    int fYSize;
    /// Minimum z value in group
    int fZMin;
    /// Maximum z value in group
    int fZMax;
    /// Difference between min and max z values in group
    int fZSize;

    /// Average position of all cells in group
    TVector3 fAveragePosition;
    /// Pad at average position of all cells in group
    trex::TTPCCellInfo3D fAveragePad;
    /// Unit volume at average position of all cells in group
    trex::TTPCUnitVolume* fAverageUnitVolume;
    /// Standard deviation of all pads in this group along x axis
    float fSigmaPadX;
    /// Standard deviation of all pads in this group along y axis
    float fSigmaPadY;
    /// Standard deviation of all pads in this group along z axis
    float fSigmaPadZ;
    /// Unique id of the pad at the average position in this group
    long fAveragePadID;
    /// Total charge contained in this group
    float fCharge;
    /// Average charge contained in this group's cells
    float fAverageCharge;

    /// Map of all hits in group
    std::map<long, trex::TTPCUnitVolume*> fHitMap;

    /// Whether internal values have yet been set
    bool fIsClosed;
    /// ID corresponding to the track end or vertex this junction is supposed to represent
    unsigned int fID;
    /// Static var holding maximum id of a group
    static unsigned int sMaxID;

    /// Calculate standard deviation of all pads in this group along all axes
    void FindSigmaPads();
    /// Calculate the unique id of the pad at the average position in this group and its position
    void FindAveragePadPosID();
    /// Calculate the total charge contained in this group
    void FindCharge();

    /// Calculate internal values
    void Close();

    /// Export class for ROOT
  //    ClassDef(TTPCVolGroup, 1);
};

#endif
