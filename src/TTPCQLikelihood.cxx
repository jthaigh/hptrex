#include "TTPCQLikelihood.hxx"
#include <TMath.h>
#include <TOARuntimeParameters.hxx>

#define PI        3.14159265358979312
#define SQRTTWOPI 2.506628274631  // TMath::Sqrt(TMath::TwoPi())
#define SQRT2     1.414213562373095 // TMath::Sqrt(2.)


TTPCQLikelihood::TTPCQLikelihood(){ 
  fMinimumPhi_Eta = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.LikFit.MinPhiEta");
}

TTPCQLikelihood::~TTPCQLikelihood(){ ;}

// *********************************************************************************
double TTPCQLikelihood::eta( double b, double phi, double sigma,double longi,double trans){
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

    if( isnan(result) ) {
	  if( ND::tpcDebug().LikFit(DB_ERROR)){
        std::cout << " Is not a Number in eta calculation in likelihood TPC recon " << std::endl; 
        std::cout << " u1 = " << u1 << "  u2 =  " << u2 << "  u3 =  " << u3 << " u4 =  " << u4 << std::endl; 
        std::cout << " l1 = " << l1 << "  l2 =  " << l2 << "  l3 =  " << l3 << " l4 =  " << l4 << std::endl; 
        std::cout << " erf(l1) = " << TMath::Erf(l1) << " erf(l2) =  " << TMath::Erf(l2) << " erf(l3) = " << TMath::Erf(l3) << " erf(l4) =  " << TMath::Erf(l4) << std::endl;
      }
    }

    result += sigma/SQRTTWOPI*(TMath::Exp(-l1*l1)-TMath::Exp(-l2*l2)
        +TMath::Exp(-l3*l3)-TMath::Exp(-l4*l4));    
  }
  else{
    double u1 = ((b+t)/SQRT2/sigma);
    double u2 = ((b-t)/SQRT2/sigma);
    result = TMath::Abs(TMath::Erf(u1)-TMath::Erf(u2));

    if( isnan(result) ) {
	  if( ND::tpcDebug().LikFit(DB_ERROR)){
        std::cout << " Is not a Number in eta calculation in likelihood TPC recon " << std::endl; 
        std::cout << " u1 = " << u1 << " u2 =   " << u2 << std::endl;
        std::cout << " erf(u1) = " << TMath::Erf(u1) << " erf(u2) =  " << TMath::Erf(u2) << std::endl;
        std::cout << " b " << b << " t " << t << " sigma " << sigma << std::endl; 
      }
    }
  }

  return result;
}     

// *********************************************************************************
// expected charge on a certain pad
double TTPCQLikelihood::q_exp( double b, double phi, double sigma,double longi,double trans){
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
