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
      void SetCell(int x, int y, int z, long id);

      /// Add a charge to the total charge of this cell
      void AddCharge(double q);
      /// Associate another hit with this cell
      void AddHit(trex::TTPCHitPad* hit);
      /// Associate multiple hits with this event
      void AddHits(std::vector< trex::TTPCHitPad* > hits);
      /// Add an event to this cell, charge weighted averaging position and incrementing charge
      void AddEvent(trex::TTPCHitPad* hit);
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

      /// Get negative peak before this hit
      double GetNegativePeakEarly(){ return fNegativePeakEarly; }
      /// Get late peak before this hit
      double GetNegativePeakLate(){ return fNegativePeakLate; }

      /// Get cell unique id
      long GetID(){ return fID;}

      /// Get distance from focus for building list
      float GetFriendDist(){ return fFriendDist;}

      // setters
      /// Set total charge in this cell
      void SetQ(int q){ fQ = q;}
      /// Set local time offset for x calculations
      void SetTimeOffset(double timeOffset){ fTimeOffset = timeOffset; }
      /// Set distance from focus for building list
      void SetFriendDist(float friendDist){ fFriendDist = friendDist;}


    private:

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

      /// Negative peak before this hit
      double fNegativePeakEarly;
      /// Late peak before this hit
      double fNegativePeakLate;

      /// Cell unique id
      long fID;

      /// Distance from focus for building list
      float fFriendDist;
  };
}

#endif
