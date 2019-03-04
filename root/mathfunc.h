#ifndef MATHFUNC_H
#define MATHFUNC_H 1

#include "TMath.h"
#include "TString.h"
#include "Riostream.h"


inline Double_t lorentzian(Double_t *x, Double_t area, Double_t mean, Double_t gamma)
{
  return area * (TMath::BreitWigner(x[0],mean,gamma));
}


inline Double_t voigth(Double_t *x, Double_t height, Double_t mean, Double_t sigma, Double_t lg)
{
  Int_t r=5;
  // To calculate the Faddeeva function with relative error less than 10^(-R).
  // R can be set by the the user subject to the constraints 2 <= R <= 5.
  return (TMath::Voigt((x[0]-mean), sigma, lg, r)) * 2.506628 * sigma * height;
}

inline Double_t voigt(Double_t *x, Double_t area, Double_t mean, Double_t sigma, Double_t lg)
{
  Int_t r=5;
  return area * (TMath::Voigt((x[0]-mean), sigma, lg, r));
}

inline Double_t letail(Double_t *x, Double_t gain, Double_t mean, Double_t sigma, Double_t beta)
// --- low energy tail distribution ---
// : the convolution of an exponential function with a Gaussian
// number of parameter ... 4
// Gain
// (Gaussian) mean : the centroid channel corresponding to the incident photon energy
// (Gaussian) sigma : the standard deviation of the Gaussian component
// Beta : the slope of the exponential feature
{
  if (gain == 0) return 0;
  if (sigma == 0) return 1.e30;
  if (beta == 0) return 1.e30;
  Double_t sqrt2 = 1.414213562373095;
  Double_t arg = (x[0] - mean)/sigma/sqrt2;
  Double_t argB1 = (x[0] - mean)/beta;
  Double_t argB2 = sigma/beta/sqrt2;
  Double_t argB3 = sigma*sigma/2./beta/beta;
  Double_t norm = 1./2./beta*exp(argB3);
  return gain * norm * exp(argB1) * TMath::Erfc(arg + argB2);
}


inline Double_t hetail(Double_t *x, Double_t gain, Double_t mean, Double_t sigma, Double_t beta)
// --- high energy tail distribution ---
// the convolution of an exponential function with a Gaussian
// # of parameter ... 4
// Gain
// (Gaussian) mean : corresponding to the incident photon energy
// (Gaussian) sigma: standard deviation of the Gaussian component
// Beta: the slope of the exponential feature
{
  if (gain == 0) return 0;
  if (sigma == 0) return 1.e30;
  if (beta == 0) return 1.e30;
  Double_t sqrt2 = 1.414213562373095;
  Double_t arg = -(x[0] - mean)/sigma/sqrt2;
  Double_t argB1 = -(x[0] - mean)/beta;
  Double_t argB2 = sigma/beta/sqrt2;
  Double_t argB3 = sigma*sigma/2./beta/beta;
  Double_t norm = 1./2./beta*exp(argB3);
  return gain * norm * exp(argB1) * TMath::Erfc(arg + argB2);
}
  


inline Double_t GetError( Double_t x1, Double_t x1err,
                          Double_t x2, Double_t x2err, TString ope="ADD")
{
  Double_t r1 = x1err / x1;
  Double_t r2 = x2err / x2;
  Double_t err=0;
  if ( (ope == "ADD") || (ope == "SUB")) {
    err = sqrt( x1err*x1err + x2err*x2err );
  } else if ( ope == "MUL" ) {
    err = (x1*x2)*sqrt( r1*r1 + r2*r2 );
  } else if ( ope == "DIV" ) {
    err = (x1/x2)*sqrt( r1*r1 + r2*r2 );
  } else {
    cout << "--- Error : 'ope' should be 'ADD','SUB','MUL','DIV' ---" << endl;
  }
  return err;
}



Double_t GetE2Ca(Double_t ene1, Double_t ene2, Double_t ch1, Double_t ch2)
{
  Double_t e2ca;
  if (ene1!=ene2) { 
    e2ca = TMath::Abs(ch2-ch1)/TMath::Abs(ene2-ene1);   
  } else {
    cout << "Error: cannot define the slope at GetE2Ca" << endl;
    exit(1);
  }
  return e2ca;
}



Double_t GetE2Cb(Double_t ene1, Double_t ene2, Double_t ch1, Double_t ch2)
{
  Double_t e2ca = GetE2Ca(ene1,ene2,ch1,ch2);
  Double_t e2cb = ch1 - e2ca*ene1;
  return e2cb;
}



Double_t GetV2Ea(Double_t ene1, Double_t ene2, Double_t vol1, Double_t vol2)
{
  Double_t v2ea;
  if (vol1!=vol2) { 
    v2ea = TMath::Abs(ene2-ene1)/TMath::Abs(vol2-vol1);   
  } else {
    cout << "Error: cannot define the slope at GetV2Ea" << endl;
    exit(1);
  }
  return v2ea;
}



Double_t GetV2Eb(Double_t ene1, Double_t ene2, Double_t vol1, Double_t vol2)
{
  Double_t v2ea = GetV2Ea(ene1,ene2,vol1,vol2);
  Double_t v2eb = ene1 - v2ea*vol1;
  return v2eb;
}


Double_t ConvV2E(Double_t vol, Double_t ene1, Double_t ene2, Double_t vol1, Double_t vol2)
{
  Double_t a = GetV2Ea(ene1,ene2,vol1,vol2);
  Double_t b = GetV2Eb(ene1,ene2,vol1,vol2);
  Double_t ene = a*vol + b;
  return ene;
}



// ---------------------------------------------------------


Double_t ConvEneFrom2(Double_t a,
		      Double_t b1, Double_t x1,
		      Double_t b2, Double_t x2,
		      Double_t x0)
{
  Double_t y1 = a*(x0-x1) + b1;
  Double_t y2 = a*(x0-x2) + b2;
  Double_t y = (y1 + y2)/2.;
  return y;
}




Double_t ConvEneErrFrom2(Double_t a, Double_t aErr,
			 Double_t b1, Double_t b1Err,
			 Double_t x1, Double_t x1Err,
			 Double_t b2, Double_t b2Err,
			 Double_t x2, Double_t x2Err,
			 Double_t x0, Double_t x0Err)
{
  // y1 = a*(x0-x1) + b1;
  // y2 = a*(x0-x2) + b2;
  
  Double_t delta_x1 = x0-x1;
  Double_t delta_x1Err = GetError(x0,x0Err,x1,x1Err,"SUB");
  Double_t mul1 = a*delta_x1;
  Double_t mul1Err = GetError(a,aErr,delta_x1,delta_x1Err,"MUL");
  Double_t y1Err = GetError(mul1,mul1Err,b1,b1Err,"ADD");

  Double_t delta_x2 = x0-x2;
  Double_t delta_x2Err = GetError(x0,x0Err,x2,x2Err,"SUB");
  Double_t mul2 = a*delta_x2;
  Double_t mul2Err = GetError(a,aErr,delta_x2,delta_x2Err,"MUL");
  Double_t y2Err = GetError(mul2,mul2Err,b2,b2Err,"ADD");

  Double_t yErr = TMath::Sqrt(y1Err*y1Err + y2Err*y2Err)/2.;

  return yErr;
}



const Double_t AvEneSilicon = 3.81;// eV @ 77K

inline Double_t GetSigmaEne(Double_t ene, Double_t noise, Double_t fano, Double_t W)
{
  // Sigma = sqrt( ConstNoise*ConstNoise + FANO * W * Ene)
  return sqrt(noise*noise + fano*W*ene);
}


inline Double_t GetSigmaCh(Double_t ch, Double_t noise, Double_t fano, Double_t W)
{
  // Sigma = sqrt( ConstNoise*ConstNoise + FANO * W * ch)
  return sqrt(noise*noise + fano*W*ch);
}


inline Double_t GetCalibSlope(Double_t Ch1, Double_t RefEne1, Double_t Ch2, Double_t RefEne2)
{
  return TMath::Abs(Ch2 - Ch1)/TMath::Abs(RefEne2 - RefEne1);
}




#endif
