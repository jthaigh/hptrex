#ifndef TTPCFEATUREFINDER3D_HXX
#define TTPCFEATUREFINDER3D_HXX

// c++
#include <iostream>
#include <cmath>
#include <map>
#include <algorithm>

// eddy
#include "TTPCLayout.hxx"

namespace ND{
  /// Class for identifying key features in a 3D array of charged cells.
  class TTPCFeatureFinder3D{
    public:
      /// Default constructor
      TTPCFeatureFinder3D(ND::TTPCLayout* layout);
      /// Default destructor
      ~TTPCFeatureFinder3D();

      /// Add charge to a specified cell
      void AddHit(int x, int y, int z, float val);
      /// Run the feature finding algorithm
      void Process();

      /// Find the cell ranges in x, y and z 
      void GetSizes(int& xSize, int& ySize, int& zSize);

      /// Return 2D array projection of fAlgMap, with input axis reduced to its maximum for each element
      float** GetHitArray(int axis=1);
      /// Create 2D array projection of fAlgMap, with input axis reduced to its maximum for each element
      float** MakeHitArray(int axis=1);

    private:
      /// Layout information for the input cells
      ND::TTPCLayout* fLayout;
      /// Layout information used during processing of input cells
      ND::TTPCLayout* fBigLayout;

      /// Map of all cells touched by algorithm
      std::map<long, float> fAlgMap;

      /// Switches whether the algorithm converts the hit list back to its original co-ordinate system or keeps expanded co-ordinates
      bool fRebaseHits;

      /// Standard deviation for kernel used to smooth input
      float fSigmaSmooth;
      /// Standard deviation over which to sample gradient
      float fSigmaStruct;
      /// Size of input area window to process in each part of smoothing
      int fWindowSmooth;
      /// Size of input area window to process in each part of gradient sampling
      int fWindowStruct;
      /// Size of input area window to process for median filter
      int fWindowMedian;
      /// Size of input area window to process for maximum filters
      int fWindowMax;

      /// Representation of z-y projection of fAlgMap
      float** fHitArrayZY;
      /// Representation of z-x projection of fAlgMap
      float** fHitArrayZX;
      /// Representation of y-x projection of fAlgMap
      float** fHitArrayYX;

      /// Smooth out input hits for feature finding
      void SmoothHits(float sigma=2.0, int window=2);
      /// Get response function for feature finding
      void GetResponse(float sigma=2.0, int window=2);

      /// Response function for feature finding
      float ResponseFunc(float i_xx, float i_yy, float i_zz, float i_xy, float i_xz, float i_yz);
      /// Original esponse function for feature finding
      float ResponseFuncClassic(float i_xx, float i_yy, float i_zz, float i_xy, float i_xz, float i_yz);
      /// Alternate response function for feature finding
      float HSResponseFunc(float i_xx, float i_yy, float i_zz, float i_xy, float i_xz, float i_yz, float k=.05);

      /// Gaussian blur or derivative filter a selection of hits, given a specified window size and standard deviation and axis
      std::map<long, float> GetFilteredHits1D(std::map<long, float> inputMap, int axis=0, int deriv=0, float sigma=2.0, int window=2, bool destructive=false);

      /// Create a gaussian kernel in 1D
      float* GaussianKernel1D(float sigma=2.0, int window=2);
      /// Create a gaussian first derivative kernel in 1D
      float* GaussianDerivativeKernel1D(float sigma=2.0, int window=2);
      /// Fast function for extracting eigenvalues from a 3x3 matrix.  From Joachim Kopp, Numerical diagnolization of hermitial 3x3 matrices, arxiv.org preprint: physics/0610206, Int. J. Mod. Phys. C19 (2008) 523-548
      int dsyevc3(double A[3][3], double w[3]);

      /// Clear all elements in a 2D array
      void ClearArray(float** arr, int sizeX);
  };
}

#endif
