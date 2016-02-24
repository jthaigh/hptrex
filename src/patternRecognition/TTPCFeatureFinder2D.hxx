#ifndef TTPCFEATUREFINDER2D_HXX
#define TTPCFEATUREFINDER2D_HXX

// c++
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

// eddy
#include "TTPCLayout.hxx"

namespace ND{
  /// Class for identifying key features in a 2D array of charged cells.
  class TTPCFeatureFinder2D{
    public:
      /// Default constructor
      TTPCFeatureFinder2D(ND::TTPCLayout* layout, int axis=1);
      /// Default destructor
      ~TTPCFeatureFinder2D();

      /// Set the minimum and maximum bounds of the cells in which to identify features directly 
      void SetRanges(int rangeX,int minX,int maxX, int rangeY,int minY,int maxY);
      /// Set up initial values and objects for feature finding
      void Initialise();
      /// Add charge to a specified cell
      void AddHit(int x, int y, float val);
      /// Run the feature finding algorithm
      void Process();

      /// Finds points at which the feature finding response function is peaked (i.e. points of interest)
      std::vector<ND::TTPCCell2D> GetLocalMaxima(int window = 3);

      /// Return array representing response function (for debugging)
      float** GetHitArray();
      /// Find the cell sizes in x and y 
      void GetSizes(int& xSize, int& ySize);

    private:
      /// Switches whether the algorithm converts the hit list back to its original co-ordinate system or keeps expanded co-ordinates
      bool fRebaseHits;
      /// Whether or not the layout has been set for this algorithm
      bool fHasLayout;

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
      /// Maximum size of input area window for any filter
      int fWindowMax;

      /// Threshold below which to discard maxima
      float fMaxThreshold;

      /// Default extent of array of input cells in x
      int fBaseSizeX;
      /// Default extent of array of input cells in y
      int fBaseSizeY;

      /// Lowest x edge of input cells
      int fMinX;
      /// Lowest y edge of input cells
      int fMinY;

      /// Highest x edge of input cells
      int fMaxX;
      /// Highest y edge of input cells
      int fMaxY;

      /// Used array of input cells in x
      int fSizeX;
      /// Used array of input cells in y
      int fSizeY;

      /// Offset to apply when feeding in input cells, in x
      int fOffsetX;
      /// Offset to apply when feeding in input cells, in y
      int fOffsetY;

      /// 2D array of cells used in all stages of processing feature finding
      float** fHitList;

      /// Gaussian smoothing of fHitList
      void SmoothHits(float sigma=2.0, int window=8);
      /// Median filter of fHitList
      void MedianFilter(int window = 3);
      /// Maximum filter of fHitList
      void MaximumFilter(int window = 3);
      /// Remove everything but local maxima from fHitList
      void LocalMaximumFilter(int window = 3);
      /// Response function on fHitList
      void GetResponse(float sigma=2.0, int window=8);
      /// Response function tuned to edge detection on fHitList
      void GetEdgeResponse(float sigma=2.0, int window=8);

      /// Return fHitList to original, input co-ordinate system
      void RebaseHits();

      /// Response function for feature finding
      float ResponseFunc(float i_xx, float i_yy, float i_xy);
      /// Alternate response function for feature finding
      float HSResponseFunc(float i_xx, float i_yy, float i_xy, float k=.05);

      /// Create a gaussian kernel
      float** GenerateKernel(float sigma=2.0, int window=8);
      /// Create a gaussian first derivative kernel along a specified axis
      float** GaussianDerivativeFilter(float sigma=2.0, int window=8, int axis=0);
      /// Create a gaussian second derivative kernel along a specified axis
      float** GaussianSecondDerivativeFilter(float sigma=2.0, int window=8, int axis=0);
      /// Create a laplacian kernel
      float** LaplacianFilter(float sigma=2.0, int window=8);
      /// Fully delete a given two dimensional array
      void ClearArray(float** arr, int sizeX);
  };
}

#endif
