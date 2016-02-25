// eddy
#include "TTPCFeatureFinder2D.hxx"

ND::TTPCFeatureFinder2D::TTPCFeatureFinder2D(ND::TTPCLayout* layout, int axis){
  fRebaseHits = true;

  // smooth across 2 pixels in both directions
  fSigmaSmooth = 3.5;
  fSigmaStruct = 3.5;
  // restrict window to the 3 sigma mark of the smoothing kernel
  fWindowSmooth = int(fSigmaSmooth * 3);
  fWindowStruct = int(fSigmaStruct * 3);
  fWindowMedian = 1;

  fMaxThreshold = 0.000001;

  // grab maxima, minima and extents in x and y from input layout object
  layout->GetRanges(fBaseSizeX,fMinX,fMaxX, fBaseSizeY,fMinY,fMaxY, axis);

  // determine maximum size of input area window
  if(fWindowSmooth < fWindowStruct) fWindowMax = fWindowStruct;
  if(fWindowMax < fWindowMedian) fWindowMax = fWindowMedian;

  // actual size in x and y is base size plus window extent above and below
  fSizeX = (fBaseSizeX) + (fWindowMax*2);
  fSizeY = (fBaseSizeY) + (fWindowMax*2);

  // offset for input plus window extent below is minimum of input
  fOffsetX = fMinX - fWindowMax;
  fOffsetY = fMinY - fWindowMax;

  // set array sizes in x and y
  const int sizeX = fSizeX;
  const int sizeY = fSizeY;

  // create 2D array of 0s
  fHitList = new float* [sizeX];
  for(int i=0; i<fSizeX; i++){
    fHitList[i] = new float [sizeY];
    for(int j=0; j<fSizeY; j++){
      fHitList[i][j] = 0.;
    };
  };
}

ND::TTPCFeatureFinder2D::~TTPCFeatureFinder2D(){
  // find size of arrays we need to delete
  int clearSize = fSizeX;
  if(fRebaseHits) clearSize = fBaseSizeX;

  // delete fHitList fully
  ClearArray(fHitList, clearSize);
}

void ND::TTPCFeatureFinder2D::SetRanges(int rangeX,int minX,int maxX, int rangeY,int minY,int maxY){
  // set maxima, minima and extents in x and y directly
  fBaseSizeX = rangeX;
  fMinX = minX;
  fMaxX = maxX;

  fBaseSizeY = rangeY;
  fMinY = minY;
  fMaxY = maxY;
}
void ND::TTPCFeatureFinder2D::AddHit(int x, int y, float val){
  // make sure nothing has gone wrong
  if (x < fMinX || y < fMinY) return;
  if (x > fMaxX || y > fMaxY) return;
  // subtract offset from input x and y and insert its value into fHitList
  fHitList[x - fOffsetX][y - fOffsetY] = val;
}
void ND::TTPCFeatureFinder2D::Process(){
  // gaussian smooth of fHitList
  SmoothHits(fSigmaSmooth, fWindowSmooth);
  // response function of fHitlist
  //GetResponse(fSigmaStruct, fWindowStruct);
  // return fHitList to input co-ordinates
  RebaseHits();
}

float** ND::TTPCFeatureFinder2D::GetHitArray(){
  return fHitList;
}
void ND::TTPCFeatureFinder2D::GetSizes(int& xSize, int& ySize){
  xSize = fBaseSizeX;
  ySize = fBaseSizeY;
}
float ND::TTPCFeatureFinder2D::ResponseFunc(float i_xx, float i_yy, float i_xy){
  // response function based in determinant and trace of structure tensor [ (i_xx) (i_xy) ][ (i_xy) (i_yy) ]
  float det = i_xx*i_yy - i_xy*i_xy;
  float tr = i_xx + i_yy;
  if(tr > 0) return det / tr;
  return 0;
}
float ND::TTPCFeatureFinder2D::HSResponseFunc(float i_xx, float i_yy, float i_xy, float k){
  // alternate response function based in determinant and trace of structure tensor [ (i_xx) (i_xy) ][ (i_xy) (i_yy) ]
  float det = i_xx*i_yy - i_xy*i_xy;
  float tr = i_xx + i_yy;
  return det - (k*tr*tr);
}
void ND::TTPCFeatureFinder2D::MedianFilter(int window){
  const int sizeX = fSizeX;
  const int sizeY = fSizeY;

  // save pointer to old hits
  float** oldHits = fHitList;

  // create and fill new 2D array
  fHitList = new float* [sizeX];
  for(int i=0; i<fSizeX; i++){
    fHitList[i] = new float [sizeY];
    for(int j=0; j<fSizeY; j++){
      // create and fill vector of every hit in window
      std::vector<int> inWindow;
      for(int i2=-window; i2<=window; i2++){
        int len = int(sqrt(window*window - i2*i2));
        for(int j2=-len; j2<=len; j2++){
          int i3 = i+i2;
          int j3 = j+j2;
          if(i3 < 0 || i3 >= sizeX || j3 < 0 || j3 >= sizeY) continue;
          inWindow.push_back(oldHits[i3][j3]);
        };
      };
      // sort vector of every hit in window
      std::sort(inWindow.begin(), inWindow.end());
      // find middle element (list should be odd size by construction)
      int entries = inWindow.size();
      int entry = entries/2;
      fHitList[i][j] = inWindow.at(entry);
    };
  };
  // delete old hits
  ClearArray(oldHits, fSizeX);
}
void ND::TTPCFeatureFinder2D::MaximumFilter(int window){
  const int sizeX = fSizeX;
  const int sizeY = fSizeY;

  // save pointer to old hits
  float** oldHits = fHitList;

  // create and fill new 2D array
  fHitList = new float* [sizeX];
  for(int i=0; i<fSizeX; i++){
    fHitList[i] = new float [sizeY];
    for(int j=0; j<fSizeY; j++){
      // find maximum value in window
      float maxVal = 0.;
      for(int i2=-window; i2<=window; i2++){
        int len = int(sqrt(window*window - i2*i2));
        for(int j2=-len; j2<=len; j2++){
          int i3 = i+i2;
          int j3 = j+j2;
          if(i3 < 0 || i3 >= sizeX || j3 < 0 || j3 >= sizeY) continue;
          maxVal = std::max(maxVal, oldHits[i3][j3]);
        };
      };
      // set this cell to maximum value in window
      fHitList[i][j] = maxVal;
    };
  };
  // delete old hits
  ClearArray(oldHits, fSizeX);
}
void ND::TTPCFeatureFinder2D::LocalMaximumFilter(int window){
  const int sizeX = fSizeX;
  const int sizeY = fSizeY;

  // save pointer to old hits
  float** oldHits = fHitList;

  // create and fill new 2D array
  fHitList = new float* [sizeX];
  for(int i=0; i<fSizeX; i++){
    fHitList[i] = new float [sizeY];
    for(int j=0; j<fSizeY; j++){
      // find maximum value in window
      float maxVal = 0.;
      for(int i2=-window; i2<=window; i2++){
        int len = int(sqrt(window*window - i2*i2));
        for(int j2=-len; j2<=len; j2++){
          int i3 = i+i2;
          int j3 = j+j2;
          if(i3 < 0 || i3 >= sizeX || j3 < 0 || j3 >= sizeY) continue;
          maxVal = std::max(maxVal, oldHits[i3][j3]);
        };
      };
      // set cell to zero if it isn't the max, its original value otherwise
      if(oldHits[i][j] == maxVal) fHitList[i][j] = maxVal;
      else fHitList[i][j] = 0.;
    };
  };
  // delete old hits
  ClearArray(oldHits, fSizeX);
}
void ND::TTPCFeatureFinder2D::SmoothHits(float sigma, int window){
  const int sizeX = fSizeX;
  const int sizeY = fSizeY;

  // create a 2D gaussian kernel
  float** kernel = GenerateKernel(sigma, window);
  // save pointer to old hits
  float** oldHits = fHitList;

  // create and fill new 2D array
  fHitList = new float* [sizeX];
  for(int i=0; i<fSizeX; i++){
    fHitList[i] = new float [sizeY];
    for(int j=0; j<fSizeY; j++){
      // convolve old hits with kernel to get new hits
      fHitList[i][j] = 0.;
      for(int i2=-window; i2<=window; i2++){
        int len = int(sqrt(window*window - i2*i2));
        for(int j2=-len; j2<=len; j2++){
          int i3 = i+i2;
          int j3 = j+j2;
          int i4 = window+i2;
          int j4 = window+j2;
          if(i3 < 0 || i3 >= sizeX || j3 < 0 || j3 >= sizeY) continue;
          fHitList[i][j] += oldHits[i3][j3] * kernel[i4][j4];
        };
      };
    };
  };
  // delete kernel
  ClearArray(kernel, window*2 + 1);
  // delete old hits
  ClearArray(oldHits, fSizeX);
}
void ND::TTPCFeatureFinder2D::GetResponse(float sigma, int window){
  const int sizeX = fSizeX;
  const int sizeY = fSizeY;

  // create a 2D gaussian derivative kernels
  float** kernelX = GaussianDerivativeFilter(sigma, window, 1);
  float** kernelY = GaussianDerivativeFilter(sigma, window, 2);
  // save pointer to old hits
  float** oldHits = fHitList;

  // create and fill new 2D array
  fHitList = new float* [sizeX];
  for(int i=0; i<fSizeX; i++){
    fHitList[i] = new float [sizeY];
    for(int j=0; j<fSizeY; j++){
      // convolve old hits with gaussian derivative kernes to get elements of structure tensor
      float i_xx = 0.;
      float i_yy = 0.;
      float i_xy = 0.;
      for(int i2=-window; i2<=window; i2++){
        int len = int(sqrt(window*window - i2*i2));
        for(int j2=-len; j2<=len; j2++){
          int i3 = i+i2;
          int j3 = j+j2;
          int i4 = window+i2;
          int j4 = window+j2;
          if(i3 < 0 || i3 >= sizeX || j3 < 0 || j3 >= sizeY) continue;
          float i_dx = oldHits[i3][j3] * kernelX[i4][j4];
          float i_dy = oldHits[i3][j3] * kernelY[i4][j4];

          i_xx += i_dx*i_dx; 
          i_yy += i_dy*i_dy;
          i_xy += i_dx*i_dy;
        };
      };
      // set elements of new hits to response function of structure tensor hits on this cell
      fHitList[i][j] = ResponseFunc(i_xx, i_yy, i_xy);
    };
  };
  // delete kernels
  ClearArray(kernelX, window*2 + 1);
  ClearArray(kernelY, window*2 + 1);
  // delete old hits
  ClearArray(oldHits, fSizeX);
}
void ND::TTPCFeatureFinder2D::GetEdgeResponse(float sigma, int window){
  const int sizeX = fSizeX;
  const int sizeY = fSizeY;

  // create a 2D laplacian kernel
  float** kernelL= LaplacianFilter(sigma, window);
  // save pointer to old hits
  float** oldHits = fHitList;

  // create and fill new 2D array
  fHitList = new float* [sizeX];
  for(int i=0; i<fSizeX; i++){
    fHitList[i] = new float [sizeY];
    for(int j=0; j<fSizeY; j++){
      // convolve old hits with laplacian kernel to get new hits
      fHitList[i][j] = 0.;
      for(int i2=-window; i2<=window; i2++){
        int len = int(sqrt(window*window - i2*i2));
        for(int j2=-len; j2<=len; j2++){
          int i3 = i+i2;
          int j3 = j+j2;
          int i4 = window+i2;
          int j4 = window+j2;
          if(i3 < 0 || i3 >= sizeX || j3 < 0 || j3 >= sizeY) continue;
          fHitList[i][j] += oldHits[i3][j3] * kernelL[i4][j4];
        };
      };
    };
  };
  // delete kernel
  ClearArray(kernelL, window*2 + 1);
  // delete old hits
  ClearArray(oldHits, fSizeX);
}
std::vector<ND::TTPCCell2D> ND::TTPCFeatureFinder2D::GetLocalMaxima(int window){
  // initialise array of maxima
  std::vector<ND::TTPCCell2D> localMaxima; 
  const int sizeX = fBaseSizeX;
  const int sizeY = fBaseSizeY;

  for(int i=0; i<sizeX; i++)
    for(int j=0; j<sizeY; j++){
      // find maximum in window (to minimum set by fMaxThreshold)
      float maxVal = fMaxThreshold;
      for(int i2=-window; i2<=window; i2++){
        int len = int(sqrt(window*window - i2*i2));
        for(int j2=-len; j2<=len; j2++){
          int i3 = i+i2;
          int j3 = j+j2;
          if(i3 < 0 || i3 >= sizeX || j3 < 0 || j3 >= sizeY) continue;
          maxVal = std::max(maxVal, fHitList[i3][j3]);
        };
      };
      // if current value is local maximum, add it to list of maxima 
      if(fHitList[i][j] == maxVal){
        ND::TTPCCell2D cell;
        cell.x = i;
        cell.y = j;
        localMaxima.push_back(cell);
      };
    };

  return localMaxima;
}

void ND::TTPCFeatureFinder2D::RebaseHits(){
  // takes arrays from large size (hit + window) to base size
  if (!fRebaseHits) return;

  const int sizeX = fBaseSizeX;
  const int sizeY = fBaseSizeY;

  // save pointer to old hits
  float** oldHits = fHitList;

  // create and fill new 2D array with same size as original input
  fHitList = new float* [sizeX];
  for(int i=0; i<sizeX; i++){
    fHitList[i] = new float [sizeY];
    for(int j=0; j<sizeY; j++){
      fHitList[i][j] = oldHits[i+fWindowMax][j+fWindowMax];
    };
  };
  // delete old hits
  ClearArray(oldHits, fSizeX);
}

float** ND::TTPCFeatureFinder2D::GenerateKernel(float sigma, int window){
  // work out size of kernel
  const int size = window*2 + 1;

  float sum=0.;

  float** kernel = new float* [size];
  kernel[window] = new float [size];
  for(int i=0; i<=window; i++){
    // set elements based from 2D gaussian function around center
    kernel[window + i] = new float [size];
    if(i > 0) kernel[window - i] = new float [size];
    for(int j=0; j<=window; j++){
      float r2 = (i*i) + (j*j);
      float val = exp(-r2/(2*sigma*sigma));
      kernel[window + i][window + j] = val;
      if(j > 0) kernel[window + i][window - j] = val;
      if(i > 0) kernel[window - i][window + j] = val;
      if(i > 0 && j > 0) kernel[window - i][window - j] = val;
      sum += val * 4.;
    };
  };
  // normalise kernel
  for(int i=-window; i<=window; i++){
    for(int j=-window; j<=window; j++){
      kernel[i+window][j+window] /= sum;
    };
  };
  return kernel;
}
float** ND::TTPCFeatureFinder2D::GaussianDerivativeFilter(float sigma, int window, int axis){
  // start from gaussian kernel
  float ** kernel = GenerateKernel(sigma, window);

  // form derivative of gaussian kernel along x or y depending on specified axis
  if(axis == 1){
    for(int i=-window; i<=window; i++){
      for(int j=-window; j<=window; j++){
        float fact = i / sigma;
        kernel[i+window][j+window] *= fact;
      };
    };
  }
  else if(axis == 2){
    for(int i=-window; i<=window; i++){
      for(int j=-window; j<=window; j++){
        float fact = j / sigma;
        kernel[i+window][j+window] *= fact;
      };
    };
  };
  return kernel;
}
float** ND::TTPCFeatureFinder2D::GaussianSecondDerivativeFilter(float sigma, int window, int axis){
  // start from gaussian kernel
  float ** kernel = GenerateKernel(sigma, window);

  // form second derivative of gaussian kernel along x or y depending on specified axis
  if(axis == 1){
    for(int i=-window; i<=window; i++){
      for(int j=-window; j<=window; j++){
        float fact = (float(i*i)/sigma - 1.) / sigma;
        kernel[i+window][j+window] *= fact;
      };
    };
  }
  else if(axis == 2){
    for(int i=-window; i<=window; i++){
      for(int j=-window; j<=window; j++){
        float fact = (float(j*j)/sigma - 1.) / sigma;
        kernel[i+window][j+window] *= fact;
      };
    };
  };

  return kernel;
}
float** ND::TTPCFeatureFinder2D::LaplacianFilter(float sigma, int window){
  // start from gaussian kernel
  float ** kernel = GenerateKernel(sigma, window);

  // form laplacian of gaussian kernel 
  for(int i=-window; i<=window; i++){
    for(int j=-window; j<=window; j++){
      float r2 = i*i + j*j;
      float fact = (r2/(sigma*sigma) - 1.) / (sigma*sigma);
      kernel[i+window][j+window] *= fact;
    };
  };
  return kernel;
}
void ND::TTPCFeatureFinder2D::ClearArray(float** arr, int sizeX){
  // delete all elements AND sub-elements of 2D array
  for(int i=0; i<sizeX; i++) delete[] arr[i];
  delete[] arr;
}
