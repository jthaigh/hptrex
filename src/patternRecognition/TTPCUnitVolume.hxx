#ifndef TTPCUNITVOLUME_HXX
#define TTPCUNITVOLUME_HXX

// c++
#include <map>
#include <vector>
#include <algorithm>

// ROOT
#include <TVector3.h>

// TREx
#include <TTPCHitPad.hxx>

// eddy
#include "TTPCLayout.hxx"

namespace trex{
  /// Class for individual volume elements used in both feature finding and path finding
  class TTPCUnitVolume{
    public:
      /// Default constructor
      TTPCUnitVolume();
      /// Default destructor
      ~TTPCUnitVolume();

      /// Set the position, edge status and unique id of this cell
      void SetCell(int x, int y, int z, int edgeX, int edgeY, int edgeZ, long id);
      /// Set other useful information on TPC location of this cell
      void SetAux(int segX, int segY, int segZ);

      /// Add a charge to the total charge of this cell
      void AddCharge(double q);
      /// Associate another hit with this cell
      void AddHit(trex::TTPCHitPad* hit);
      /// Associate multiple hits with this event
      void AddHits(std::vector< trex::TTPCHitPad* > hits);
      /// Add an event to this cell, charge weighted averaging position and incrementing charge
      void AddEvent(trex::TTPCHitPad* hit);
      /// Get edge status of this cell in x, y and z - -1 for a lower edge of a MM volume, 1 for an upper edge, -2 for a dead region
      void GetEdges(int& edgeX, int& edgeY, int& edgeZ);
      /// Summarise position and edge status in one struct
      trex::TTPCCellInfo3D GetCellInfo3D();

      /// Get total number of peaks
      int GetNPeaksSum();
      /// Get maximum number of peaks
      int GetNPeaksMax();
      /// Get number of saturated peaks
      int GetNSaturated();
      /// Get total saturation
      int GetSaturation();

      /// Reset friends for association
      void ClearExtendedCell(bool isFocus=true);

      // getters
      /// Get whether or not position is set
      bool GetHasPos(){ return fHasPos;}

      /// Get cell x id
      int GetX(){ return fX;}
      /// Get cell y id
      int GetY(){ return fY;}
      /// Get cell z id
      int GetZ(){ return fZ;}

      /// Get total charge in this cell
      double GetQ(){ return fQ;}
      /// Get maxium charge in this cell
      double GetQMax(){ return fQMax;}
      /// Get average position in this cell
      TVector3 GetPos(){ return fPos;}
      /// Get cell index representation ov average position in this cell
      TVector3 GetPosXYZ(){ return TVector3( (double)fX, (double)fY, (double)fZ); }
      /// Get average time in this cell
      double GetTime(){ return fTime;}
      /// Get offset average time in this cell
      double GetTimeNom(){ return fTimeNom;}
      /// Get minimum time in this cell
      double GetTimeMin(){ return fTimeMin;}
      /// Get maximum time in this cell
      double GetTimeMax(){ return fTimeMax;}

      /// Get hits associated with this cell
      std::vector< trex::TTPCHitPad* > GetHits(){ return fHits;}
      /// Get iterator to beginning of hits associated with this cell
      std::vector< trex::TTPCHitPad* >::iterator GetHitsBegin(){ return fHits.begin();}
      /// Get iterator to end of hits associated with this cell
      std::vector< trex::TTPCHitPad* >::iterator GetHitsEnd(){ return fHits.end();}
      /// Get size of list of hits associated with this cell
      unsigned int size(){ return fHits.size();}

      /// Get cell x edge status
      int GetEdgeX(){ return fEdgeX;}
      /// Get cell y edge status
      int GetEdgeY(){ return fEdgeY;}
      /// Get cell z edge status
      int GetEdgeZ(){ return fEdgeZ;}

      /// Get TPC half the cell is in
      int GetSegX(){ return fSegX; }
      /// Get MM row ID for the cell
      int GetSegY(){ return fSegY; }
      /// Get MM col ID for the cell
      int GetSegZ(){ return fSegZ; }

      /// Get negative peak before this hit
      double GetNegativePeakEarly(){ return fNegativePeakEarly; }
      /// Get late peak before this hit
      double GetNegativePeakLate(){ return fNegativePeakLate; }

      /// Get cell unique id
      long GetID(){ return fID;}
      /// Get TPC for this cell
      unsigned int GetTPC(){ return fTPC; }
      /// Get half for this cell
      unsigned int GetHalf(){ return fHalf; }
      /// Get MicroMega for this cell
      unsigned int GetMM(){ return fMM; }

      /// Get distance from focus for building list
      float GetFriendDist(){ return fFriendDist;}

      // setters
      /// Set total charge in this cell
      void SetQ(int q){ fQ = q;}
      /// Set local time offset for x calculations
      void SetTimeOffset(double timeOffset){ fTimeOffset = timeOffset; }
      /// Set MM location
      void SetMMLoc(unsigned int tpc, unsigned int half, unsigned int mm){ fTPC = tpc; fHalf = half; fMM = mm; }
      /// Set distance from focus for building list
      void SetFriendDist(float friendDist){ fFriendDist = friendDist;}


    private:
      /// Total number of TPCs
      static const int scNTPC;
      /// Total number of halves
      static const int scNHalf;
      /// Total number of MMs
      static const int scNMM;

      /// Local time offset for x calculations
      double fTimeOffset;

      /// Whether or not position is set
      bool fHasPos;

      /// Cell x id
      int fX;
      /// Cell y id
      int fY;
      /// Cell z id
      int fZ;

      /// Total charge in this cell
      double fQ;
      /// Maximum charge in this cell
      double fQMax;
      /// Average position in this cell
      TVector3 fPos;
      /// Average time in this cell
      double fTime;
      /// Average offset time in this cell
      double fTimeNom;
      /// Minimum time in this cell
      double fTimeMin;
      /// Maximum time in this cell
      double fTimeMax;

      /// Hits associated with this cell
      std::vector< trex::TTPCHitPad* > fHits;

      /// Cell x edge status
      int fEdgeX;
      /// Cell y edge status
      int fEdgeY;
      /// Cell z edge status
      int fEdgeZ;

      /// TPC half the cell is in
      int fSegX;
      /// MM row ID for the cell
      int fSegY;
      /// MM col ID for the cell
      int fSegZ;

      /// Negative peak before this hit
      double fNegativePeakEarly;
      /// Late peak before this hit
      double fNegativePeakLate;

      /// Cell unique id
      long fID;
      /// TPC for this cell
      unsigned int fTPC;
      /// Half for this cell
      unsigned int fHalf;
      /// MicroMega for this cell
      unsigned int fMM;

      /// Distance from focus for building list
      float fFriendDist;
  };
}

#endif
