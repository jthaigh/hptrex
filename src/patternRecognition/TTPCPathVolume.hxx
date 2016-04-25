#ifndef TTPCPATHVOLUME_HXX
#define TTPCPATHVOLUME_HXX

// c++
#include <stdexcept>
#include <vector>
#include <algorithm>

// eddy
#include "TTPCUnitVolume.hxx"

namespace trex{
  // forward declaration
  class TTPCPathVolume;

  /// Class for elements of volume associated with a path
  class TTPCPathVolume{
    public:
      /// Constructor with unit volume attached
      TTPCPathVolume(trex::TTPCUnitVolume* unitVolume=0);
      /// Default destructor
      ~TTPCPathVolume();

      /// Get this and cells in its assocaited cluster.  Filter of 1 for only vertical, 2 for only horizontal.
      std::vector<trex::TTPCUnitVolume*> GetExtendedCell(int filter=0);
      /// Get hits associated with associated cell and all connected to it
      std::vector<trex::TTPCHitPad*> GetHits();

      // getters
      /// Get the unit volume associated with this path volume
      trex::TTPCUnitVolume* GetUnitVolume(){ return fUnitVolume;}
      /// Get indicator of whether list of friends linked to this contains a given volume
      bool GetFriendsContains(trex::TTPCUnitVolume* vol){ return (std::find(fFriends.begin(), fFriends.end(), vol) != fFriends.end()); }
      /// Get iterator to start of list of friends linked to this
      std::vector<trex::TTPCUnitVolume*>::iterator GetFriendsBegin(){ return fFriends.begin();}
      /// Get iterator to end of list of friends linked to this
      std::vector<trex::TTPCUnitVolume*>::iterator GetFriendsEnd(){ return fFriends.end();}

      /// Get whether this has elements attached
      bool GetHasCluster(){ return fFriends.size()>0; }
      /// Get number of elements attached
      int GetClusterSize(){ return (int)fFriends.size(); }

      /// Get whether or not associated cell's position is set
      bool GetHasPos(){ return fUnitVolume->GetHasPos();}

      /// Get associated cell x id
      int GetX(){ return fUnitVolume->GetX();}
      /// Get associated cell y id
      int GetY(){ return fUnitVolume->GetY();}
      /// Get associated cell z id
      int GetZ(){ return fUnitVolume->GetZ();}

      /// Get total charge in associated cell
      int GetQ(){ return fUnitVolume->GetQ();}
      /// Get average position in associated cell
      TVector3 GetPos(){ return fUnitVolume->GetPos();}
      /// Get average position in terms of xyz cell of associated cell
      TVector3 GetPosXYZ(){ return fUnitVolume->GetPosXYZ();}

      /// Get average position in terms of xyz of the whole cluster
      TVector3 GetAvgPosXYZ();

      /// Get hits associated with associated cell
      std::vector< trex::TTPCHitPad* > GetCellHits(){ return fUnitVolume->GetHits();}

      /// Set minimum and maximum values
      void Close();

      /// Get minimum x cell
      int GetXMin(){ Close(); return fXMin;}
      /// Get maximum x cell
      int GetXMax(){ Close(); return fXMax;}
      /// Get spread in x
      int GetXSize(){ Close(); return fXSize;}
      /// Get minimum y cell
      int GetYMin(){ Close(); return fYMin;}
      /// Get maximum y cell
      int GetYMax(){ Close(); return fYMax;}
      /// Get spread in y
      int GetYSize(){ Close(); return fYSize;}
      /// Get minimum z cell
      int GetZMin(){ Close(); return fZMin;}
      /// Get maximum z cell
      int GetZMax(){ Close(); return fZMax;}
      /// Get spread in z
      int GetZSize(){ Close(); return fZSize;}
      /// Get average position
      TVector3 GetAveragePos(){ Close(); return fAveragePos; }
      /// Get average position in x, y and z
      TVector3 GetAveragePosXYZ(){ Close(); return fAveragePosXYZ; }

      /// Get associated cell unique id
      long GetID(){ return fUnitVolume->GetID();}

      /// Get whether this is a vertical cluster
      bool GetIsVertical(){ return fIsVertical;}
      /// Get flag for whether this is an x-cluster
      bool GetIsXCluster(){ return fIsXCluster;}
      /// Get y-z angle determined by pattern recognition (currently just ranges from 0-90 degrees)
      float GetPatRecAngle(){ return fPatRecAngle;}

      /// Mark friend for deletion
      void MarkFriend(std::vector<trex::TTPCUnitVolume*>::iterator focusFriendIt);
      /// Clear friends marked for deletion
      void ClearMarked();
      /// Clear list of friends
      void ClearFriends();

      // setters
      /// Add new friend to list of friends linked to this
      void AddFriend(trex::TTPCUnitVolume* focusFriend){ fClosed=false; fFriends.push_back(focusFriend); }
      /// Add new friend to list of friends linked to this, checking to avoid duplicates
      void SafeAddFriend(trex::TTPCUnitVolume* focusFriend){ if(!GetFriendsContains(focusFriend)) AddFriend(focusFriend); }
      /// Set the unit volume associated with this path volume
      void SetUnitVolume(trex::TTPCUnitVolume* unitVolume){ fClosed=false; fUnitVolume = unitVolume; }
      /// Set flag for whether this is a vertical cluster
      void SetIsVertical(bool isVertical){ fIsVertical = isVertical; }
      /// Flag for whether this is an x-cluster
      void SetIsXCluster(bool isXCluster){ fIsXCluster = isXCluster; }
      /// Set y-z angle determined by pattern recognition (currently just ranges from 0-90 degrees)
      void SetPatRecAngle(float patRecAngle){ fPatRecAngle = patRecAngle;}

      /// Print the positions of every hit in the cluster
      void PrintPositions(bool showOrientation=false);
      /// Debugging tool for checking hits for memory issues
      void CheckHits();

    private:
      /// Unit of volume associated with this element
      trex::TTPCUnitVolume* fUnitVolume;

      /// List of members of same cluster as this
      std::vector<trex::TTPCUnitVolume*> fFriends;

      /// Whether internal values have been calculated
      bool fClosed;
      /// Minimum x cell
      int fXMin;
      /// Maximum x cell
      int fXMax;
      /// Spread in x
      int fXSize;
      /// Minimum y cell
      int fYMin;
      /// Maximum y cell
      int fYMax;
      /// Spread in y
      int fYSize;
      /// Minimum z cell
      int fZMin;
      /// Maximum z cell
      int fZMax;
      /// Spread in z
      int fZSize;
      /// Average position
      TVector3 fAveragePos;
      /// Average position in x, y and z
      TVector3 fAveragePosXYZ;

      /// Flag for whether this is a vertical cluster
      bool fIsVertical;
      /// Flag for whether this is an x-cluster
      bool fIsXCluster;
      /// y-z angle determined by pattern recognition (currently just ranges from 0-90 degrees)
      float fPatRecAngle;
  };
}

#endif
