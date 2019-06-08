#ifndef FITPARAM_HH
#define FITPARAM_HH 1

#include "TString.h"

#include <iostream>
#include <fstream>

struct fitparam
{
public:
  // be careful, these are public member parameters
  TString element,linetype;
  Bool_t LET,HET,BG;
  Int_t npar;
  Double_t lrange,hrange;
  Double_t e2c,e2cErr;// [ev/ch] derivative of spline3 function
  Double_t ampl,amplErr;// amplitude of the first peak
  Double_t mean,meanErr;// mean of the first peak
  Double_t gaussian_sigma,gaussian_sigmaErr;// energy resolution of the first peak
  Double_t letail_ratio,letail_ratioErr;// low energy tail ratio
  Double_t letail_beta,letail_betaErr;// [eV] low energy tail beta
  Double_t hetail_ratio,hetail_ratioErr;// high energy tail ratio
  Double_t hetail_beta,hetail_betaErr;// [eV] high energy tail beta
  Double_t bg_0,bg_0Err;// constant value of pol1 background
  Double_t bg_1,bg_1Err;// slope value  of pol1 background
  Double_t ene,eneErr;
  Double_t enediff,enediffErr;
  Double_t fwhm,fwhmErr;
  Double_t chisq,ndf,nchi;
  Double_t width,widthErr;// physical width (Gamma)

  fitparam(){}
  fitparam(const TString finame) {
    // read a binary file and get the structure
    fitparam fp={};// this is crazy?
    ifstream file(finame, ios::in | ios::binary);
    if (!file) {
      cerr << "Error: cannot open file " << finame << endl;
      exit(1);
    }else{
      file.read((char*)this, sizeof(fp));
      cout << "reading file " << finame << endl;
      //this->Print();
    }
  }
  ~fitparam(){}
  
  void Print(Bool_t ENE=false)
  {
    cout << endl;    
    cout << Form("### FitXray %s%s ###",element.Data(),linetype.Data()) << endl;
    cout << Form("Npar = %d",npar) << endl;
    cout << Form("Chisqure/Ndf = %9.4f / %.1f = %6.4f",chisq,ndf,nchi) << endl;
    cout << Form("Fit range = (%6.1f, %6.1f)",lrange,hrange) << endl;
    cout << " --- Voigt --- " << endl;
    cout << Form(" e2c = %9.4f +- %6.4f (eV/ch)",e2c,e2cErr) << endl;
    cout << Form("area = %.5e  +- %6.4f",ampl,amplErr) << endl;
    if (!ENE) {cout << Form("mean = %9.4f +- %6.4f (ch)",mean,meanErr) << endl;}
    else {
      cout << Form("mean = %9.4f +- %6.4f (eV)",mean,meanErr) << endl;
      cout << Form("diff = %9.4f +- %6.4f (eV)",enediff,enediffErr) << endl;
    }
    cout << Form("FWHM = %9.4f +- %6.4f (eV)",fwhm,fwhmErr) << endl;
    if (LET) {
      cout << " --- Low-energy tail --- " << endl;
      cout << Form("tail ratio = %9.5f +- %6.5f",letail_ratio,letail_ratioErr) << endl;
      cout << Form("tail beta  = %9.5f +- %6.5f",letail_beta,letail_betaErr) << endl;
    }
    if (HET) {
      cout << " --- High-energy tail --- " << endl;
      cout << Form("tail ratio = %9.5f +- %6.5f",hetail_ratio,hetail_ratioErr) << endl;
      cout << Form("tail beta  = %9.5f +- %6.5f",hetail_beta,hetail_betaErr) << endl;
    }
    if (BG) {
      cout << " --- BG --- " << endl;
      cout << Form("const = %9.5f +- %6.5f",bg_0,bg_0Err) << endl;
      cout << Form("slope = %9.5f +- %6.5f",bg_1,bg_1Err) << endl;
    }
    cout << endl;    
  }
};



inline fitparam get_fp(TString finame)
{
  // how to get fitparam from binary
  // finame: binary file name
  fitparam fp(finame);
  return fp;
}


#endif
