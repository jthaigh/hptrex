// eddy
#include "TTPCFeatureFinder3D.hxx"

ND::TTPCFeatureFinder3D::TTPCFeatureFinder3D(ND::TTPCLayout* layout){
  // initialise empty map of hits
  fAlgMap = std::map<long, float>();

  fHitArrayZY = 0;

  fHitArrayZX = 0;
  fHitArrayYX = 0;

  fSigmaSmooth = 5.;
  fSigmaStruct = 5.;

  float windowMult = 2.0;
  fWindowSmooth = int(windowMult*fSigmaSmooth);
  fWindowStruct = int(windowMult*fSigmaStruct);
  fWindowMedian = std::max(fWindowSmooth, fWindowStruct);

  // assign input layout and create new layout for feature finder of larger size to accomodate smoothing
  fLayout = layout;
  fBigLayout = new ND::TTPCLayout();

  // determine maximum of all filter windows
  fWindowMax = fWindowSmooth;
  if(fWindowSmooth < fWindowStruct) fWindowMax = fWindowStruct;
  if(fWindowMax < fWindowMedian) fWindowMax = fWindowMedian;

  // set dimensions of new layout, expanded from initial layout to accomodate smoothing
  fBigLayout->SetTimeRanges(fLayout->GetTNMin(), fLayout->GetTNMax(), fLayout->GetTPMin(), fLayout->GetTPMax());
  fBigLayout->SetRanges(fLayout->GetMinX()-fWindowMax,fLayout->GetMaxX()+fWindowMax, fLayout->GetMinY()-fWindowMax,fLayout->GetMaxY()+fWindowMax, fLayout->GetMinZ()-fWindowMax,fLayout->GetMaxZ()+fWindowMax);
}
ND::TTPCFeatureFinder3D::~TTPCFeatureFinder3D(){
  // remove feature finder layout
  delete fBigLayout;

  // delete all 2D projections
  if(fHitArrayZY) ClearArray(fHitArrayZY, fLayout->GetSizeZ());
  if(fHitArrayZX) ClearArray(fHitArrayZX, fLayout->GetSizeZ());
  if(fHitArrayYX) ClearArray(fHitArrayYX, fLayout->GetSizeY());
}

void ND::TTPCFeatureFinder3D::AddHit(int x, int y, int z, float val){
  // create box of 0 initialised cells around input to accomodate smoothing
  for(int i=x-fWindowMax; i<=x+fWindowMax; i++)
    for(int j=y-fWindowMax; j<=y+fWindowMax; j++)
      for(int k=z-fWindowMax; k<=z+fWindowMax; k++){
        long id = fBigLayout->Mash(i, j, k);

        if(fAlgMap.find(id) == fAlgMap.end()) fAlgMap[id] = 0.;
      };

  // add value of current cell to its correct position
  long id = fBigLayout->Mash(x, y, z);
  fAlgMap[id] += val;
}
void ND::TTPCFeatureFinder3D::Process(){
  SmoothHits(fSigmaSmooth, fWindowSmooth);
  GetResponse(fSigmaStruct, fWindowStruct);
}

void ND::TTPCFeatureFinder3D::SmoothHits(float sigma, int window){
  // smooth map along x, y and z
  fAlgMap = GetFilteredHits1D(fAlgMap, 0, 0, sigma, window, true);
  fAlgMap = GetFilteredHits1D(fAlgMap, 1, 0, sigma, window, true);
  fAlgMap = GetFilteredHits1D(fAlgMap, 2, 0, sigma, window, true);
}
void ND::TTPCFeatureFinder3D::GetResponse(float sigma, int window){
  // get x, y and z derivative filters
  std::map<long, float> xDeriv;
  std::map<long, float> yDeriv;
  std::map<long, float> zDeriv;

  // smooth along y and z, get derivative along x
  xDeriv = GetFilteredHits1D(fAlgMap, 1, 0, sigma, window, false);
  xDeriv = GetFilteredHits1D(xDeriv, 2, 0, sigma, window, true);
  //xDeriv = GetFilteredHits1D(xDeriv, 0, 1, sigma, window, true);

  // smooth along z and x, get derivative along y
  yDeriv = GetFilteredHits1D(fAlgMap, 2, 0, sigma, window, false);
  yDeriv = GetFilteredHits1D(yDeriv, 0, 0, sigma, window, true);
  //yDeriv = GetFilteredHits1D(yDeriv, 1, 1, sigma, window, true);

  // smooth along x and y, get derivative along z
  zDeriv = GetFilteredHits1D(fAlgMap, 0, 0, sigma, window, false);
  zDeriv = GetFilteredHits1D(zDeriv, 1, 0, sigma, window, true);
  //zDeriv = GetFilteredHits1D(zDeriv, 2, 1, sigma, window, true);

  // colvolve existing elements with gaussian
  for(std::map<long, float>::iterator el = fAlgMap.begin(); el != fAlgMap.end(); ++el){
    long id = el->first;

    float dx = xDeriv[id];
    float dy = yDeriv[id];
    float dz = zDeriv[id];

    // initialise structure tensor elements to zero 
    float i_xx = dx*dx;
    float i_yy = dy*dy;
    float i_zz = dz*dz;
    float i_xy = dx*dy;
    float i_xz = dx*dz;
    float i_yz = dy*dz;

    fAlgMap[el->first] = ResponseFuncClassic(i_xx, i_yy, i_zz, i_xy, i_xz, i_yz);
  };

  xDeriv.clear();
  yDeriv.clear();
  zDeriv.clear();
}

float ND::TTPCFeatureFinder3D::ResponseFunc(float i_xx, float i_yy, float i_zz, float i_xy, float i_xz, float i_yz){
  // response function based in determinant and trace of structure tensor
  double strucTensor[3][3] = 
  {   {i_xx, i_xy, i_xz},
    {i_xy, i_yy, i_yz},
    {i_xz, i_yz, i_zz}    };
  double eigenvalues[3];
  dsyevc3(strucTensor, eigenvalues);

  double big = std::max(eigenvalues[0], eigenvalues[1]);
  big = std::max(big, eigenvalues[2]);

  return big / (eigenvalues[0] + eigenvalues[1] + eigenvalues[2]);
}
float ND::TTPCFeatureFinder3D::ResponseFuncClassic(float i_xx, float i_yy, float i_zz, float i_xy, float i_xz, float i_yz){
  // response function based in determinant and trace of structure tensor
  float det = i_xx*(i_xx*i_zz - i_yz*i_yz) - i_xy*(i_xy*i_zz - i_yz*i_xz) + i_xz*(i_xy*i_yz - i_yy*i_xz);;
  float tr = i_xx + i_yy + i_zz;
  if(tr > 0) return det / tr;
  return 0;
}
float ND::TTPCFeatureFinder3D::HSResponseFunc(float i_xx, float i_yy, float i_zz, float i_xy, float i_xz, float i_yz, float k){
  // alternate response function based in determinant and trace of structure tensor
  float det = i_xx*(i_xx*i_zz - i_yz*i_yz) - i_xy*(i_xy*i_zz - i_yz*i_xz) + i_xz*(i_xy*i_yz - i_yy*i_xz);;
  float tr = i_xx + i_yy + i_zz;
  return det - (k*tr*tr);
}

std::map<long, float> ND::TTPCFeatureFinder3D::GetFilteredHits1D(std::map<long, float> inputMap, int axis, int deriv, float sigma, int window, bool destructive){
  float* kernel;

  if(deriv == 0) kernel = GaussianKernel1D(sigma, window);
  else kernel = GaussianDerivativeKernel1D(sigma, window);

  std::map<long, float> newMap = std::map<long, float>();

  // calculate each element from sum of neighbors 
  for(std::map<long, float>::iterator el = inputMap.begin(); el != inputMap.end(); ++el){
    ND::TTPCCell3D cell = fBigLayout->UnMash(el->first);

    // initialise new element to zero
    newMap[el->first] = 0.;
    // set up starting cells for sum
    int cellX = cell.x;
    int cellY = cell.y;
    int cellZ = cell.z;
    if(axis == 0) cellX -= window+1;
    else if(axis == 1) cellY -= window+1;
    else if(axis == 2) cellZ -= window+1;
    // sum based on kernel elements
    for(int i=-window; i<=window; i++){
      // increment cell along appropriate axis
      if(axis == 0) cellX ++;
      else if(axis == 1) cellY ++;
      else cellZ ++;

      long id = fBigLayout->SafeMash(cellX, cellY, cellZ);
      if(id < 0) continue;

      std::map<long, float>::iterator it = inputMap.find(id);
      if(it == inputMap.end()) continue;

      newMap[el->first] += inputMap[id] * kernel[i+window];
    };
  };

  delete kernel;
  if(destructive) inputMap.clear();

  return newMap;
}

float* ND::TTPCFeatureFinder3D::GaussianKernel1D(float sigma, int window){
  // work out size of kernel
  const int size = window*2 + 1;

  float sum=0.;

  float* kernel = new float [size];
  for(int i=0; i<size; i++){
    int i_2 = i-window;
    int r2 = (i_2*i_2);
    float val = exp(-r2 / (2*sigma*sigma));

    kernel[i] = val;
    sum += val;
  };


  // normalise kernel
  for(int i=0; i<size; i++)
    kernel[i] /= sum;

  return kernel;
}
float* ND::TTPCFeatureFinder3D::GaussianDerivativeKernel1D(float sigma, int window){
  // start from gaussian kernel
  float* kernel = GaussianKernel1D(sigma, window);

  // work out size of kernel
  const int size = window*2 + 1;

  // form derivative of gaussian kernel 
  for(int i=0; i<size; i++){
    float fact = (i - window) / sigma;
    kernel[i] *= fact;
  };
  return kernel;
}

void ND::TTPCFeatureFinder3D::GetSizes(int& xSize, int& ySize, int& zSize){
  xSize = fLayout->GetSizeX();
  ySize = fLayout->GetSizeY();
  zSize = fLayout->GetSizeZ();
}

float** ND::TTPCFeatureFinder3D::GetHitArray(int axis){
  // for each axis, if corresponding projection is undefined create it; in both cases then return it
  if(axis == 1){
    if(!fHitArrayZY) fHitArrayZY = MakeHitArray(1);
    return fHitArrayZY;
  }
  else if(axis == 2){
    if(!fHitArrayZX) fHitArrayZX = MakeHitArray(2);
    return fHitArrayZX;
  }
  else{
    if(!fHitArrayYX) fHitArrayYX = MakeHitArray(3);
    return fHitArrayYX;
  };
}
float** ND::TTPCFeatureFinder3D::MakeHitArray(int axis){
  int sizeX = 0;
  int sizeY = 0;

  // set dimensions of output array based on axis specified for projection
  if(axis == 1){
    sizeX = fLayout->GetSizeZ();
    sizeY = fLayout->GetSizeY();
  }
  else if(axis == 2){
    sizeX = fLayout->GetSizeZ();
    sizeY = fLayout->GetSizeX();
  }
  else{
    sizeX = fLayout->GetSizeY();
    sizeY = fLayout->GetSizeX();
  };

  const int sizeI = sizeX;
  const int sizeJ = sizeY;

  // create output array of 0s
  float** hitArray = new float* [sizeI];
  for(int i=0; i<sizeI; i++){
    hitArray[i] = new float [sizeJ];
    for(int j=0; j<sizeJ; j++){
      hitArray[i][j] = 0.;
    };
  };

  // loop over all elements in 3D space
  for(std::map<long, float>::iterator el = fAlgMap.begin(); el != fAlgMap.end(); ++el){
    ND::TTPCCell3D cell = fBigLayout->UnMash(el->first);

    int i = cell.x - fLayout->GetMinX();
    int j = cell.y - fLayout->GetMinY();
    int k = cell.z - fLayout->GetMinZ(); 

    if(i < 0 || j < 0 || k < 0) continue;
    if(i >= fLayout->GetSizeX() || j >= fLayout->GetSizeY() || k >= fLayout->GetSizeZ()) continue;

    if(axis == 1){
      hitArray[k][j] = std::max(hitArray[k][j], el->second);
    }
    else if(axis == 2){
      hitArray[k][i] = std::max(hitArray[k][i], el->second);
    }
    else {
      hitArray[j][i] = std::max(hitArray[j][i], el->second);
    };
  };

  return hitArray;
}

int ND::TTPCFeatureFinder3D::dsyevc3(double A[3][3], double w[3])
  // ----------------------------------------------------------------------------
  // Calculates the eigenvalues of a symmetric 3x3 matrix A using Cardano's
  // analytical algorithm.
  // Only the diagonal and upper triangular parts of A are accessed. The access
  // is read-only.
  // ----------------------------------------------------------------------------
  // Parameters:
  //   A: The symmetric input matrix
  //   w: Storage buffer for eigenvalues
  // ----------------------------------------------------------------------------
  // Return value:
  //   0: Success
  //  -1: Error
  // ----------------------------------------------------------------------------
{
  double m, c1, c0;

  // Determine coefficients of characteristic poynomial. We write
  //       | a   d   f  |
  //  A =  | d*  b   e  |
  //       | f*  e*  c  |
  double de = A[0][1] * A[1][2];                                    // d * e
  double dd = A[0][1]*A[0][1];                                         // d^2
  double ee = A[1][2]*A[1][2];                                         // e^2
  double ff = A[0][2]*A[0][2];                                         // f^2
  m  = A[0][0] + A[1][1] + A[2][2];
  c1 = (A[0][0]*A[1][1] + A[0][0]*A[2][2] + A[1][1]*A[2][2])        // a*b + a*c + b*c - d^2 - e^2 - f^2
    - (dd + ee + ff);
  c0 = A[2][2]*dd + A[0][0]*ee + A[1][1]*ff - A[0][0]*A[1][1]*A[2][2]
    - 2.0 * A[0][2]*de;                                     // c*d^2 + a*e^2 + b*f^2 - a*b*c - 2*f*d*e)

  double p, sqrt_p, q, c, s, phi;
  const double M_SQRT3 = 1.73205080756887729352744634151;

  p = m*m - 3.0*c1;
  q = m*(p - (3.0/2.0)*c1) - (27.0/2.0)*c0;
  sqrt_p = sqrt(fabs(p));

  phi = 27.0 * ( 0.25*c1*c1*(p - c1) + c0*(q + 27.0/4.0*c0));
  phi = (1.0/3.0) * atan2(sqrt(fabs(phi)), q);

  c = sqrt_p*cos(phi);
  s = (1.0/M_SQRT3)*sqrt_p*sin(phi);

  w[1]  = (1.0/3.0)*(m - c);
  w[2]  = w[1] + s;
  w[0]  = w[1] + c;
  w[1] -= s;

  return 0;
}

void ND::TTPCFeatureFinder3D::ClearArray(float** arr, int sizeX){
  for(int i=0; i<sizeX; i++) delete arr[i];
  delete arr;
}
