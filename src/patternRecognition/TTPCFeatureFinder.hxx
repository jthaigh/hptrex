#ifndef TTPCFEATUREFINDER_HXX
#define TTPCFEATUREFINDER_HXX

// c++
#include <map>

// eddy
#include "TTPCLayout.hxx"
#include "TTPCUnitVolume.hxx"
#include "TTPCFeatureFinder2D.hxx"
#include "TTPCFeatureFinder3D.hxx"

namespace ND{
  class TTPCFeatureFinder{
    public:
      TTPCFeatureFinder(ND::TTPCLayout* layout);
      ~TTPCFeatureFinder();

      void AddHits(std::map<long, ND::TTPCUnitVolume*>);
      void Process();

      std::vector<ND::TTPCCell3D> GetFeatures(int closeness=10, int window=10);

      float** GetHitArray(int view);
      void GetSizes(int& xSize, int& ySize, int& zSize);

    private:
      bool fUseCharge;
      bool fUse3D;

      ND::TTPCLayout* fLayout;

      int fSizeX;
      int fSizeY;
      int fSizeZ;

      int fMinX;
      int fMinY;
      int fMinZ;

      int fMaxX;
      int fMaxY;
      int fMaxZ;

      int sizeX;
      int sizeY;

      ND::TTPCFeatureFinder2D* fFeaturesZY;
      ND::TTPCFeatureFinder2D* fFeaturesZX;
      ND::TTPCFeatureFinder2D* fFeaturesYX;

      ND::TTPCFeatureFinder3D* fFeatures3D;
  };
}

#endif
