#include "TTPCQLikelihood.hxx"
#include <TMath.h>

#define PI        3.14159265358979312
#define SQRTTWOPI 2.506628274631  // TMath::Sqrt(TMath::TwoPi())
#define SQRT2     1.414213562373095 // TMath::Sqrt(2.)


trex::TTPCQLikelihood::TTPCQLikelihood(){

  fMinimumPhi_Eta = 1.e-20;
}

trex::TTPCQLikelihood::~TTPCQLikelihood(){ ;}

// *********************************************************************************
double trex::TTPCQLikelihood::eta( double b, double phi, double sigma,double longi,double trans){
  double result;
  double cosphi = TMath::Cos(phi);;
  double hsinphi = longi*0.5*TMath::Abs(TMath::Sin(phi));
  double t = trans*0.5;

  if( TMath::Abs(phi) > fMinimumPhi_Eta ) {

    double u1 = ((b+t)*cosphi+hsinphi);
    double u2 = ((b+t)*cosphi-hsinphi);
    double u3 = ((b-t)*cosphi-hsinphi);
    double u4 = ((b-t)*cosphi+hsinphi);

    double factsigmasqrt2 = (SQRT2*sigma);

    double l1 = u1/factsigmasqrt2;
    double l2 = u2/factsigmasqrt2;
    double l3 = u3/factsigmasqrt2;
    double l4 = u4/factsigmasqrt2;

    result = 0.5*(  u1*TMath::Erf(l1)- u2*TMath::Erf(l2)
        + u3*TMath::Erf(l3)- u4*TMath::Erf(l4));

    result += sigma/SQRTTWOPI*(TMath::Exp(-l1*l1)-TMath::Exp(-l2*l2)
        +TMath::Exp(-l3*l3)-TMath::Exp(-l4*l4));    
  }
  else{
    double u1 = ((b+t)/SQRT2/sigma);
    double u2 = ((b-t)/SQRT2/sigma);
    result = TMath::Abs(TMath::Erf(u1)-TMath::Erf(u2));

  }

  return result;
}     

// *********************************************************************************
// expected charge on a certain pad
double trex::TTPCQLikelihood::q_exp( double b, double phi, double sigma,double longi,double trans){
  double result;

  result = eta(b,phi,sigma,longi,trans);

  if( result < 0 ) 
  {
    if( phi < 0 ) 
      result = -result; 
    else 
      result = 0.0; 
  }

  return result;
}
