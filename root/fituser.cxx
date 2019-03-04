#ifndef FITELEM_CXX
#define FITELEM_CXX 1

#include <vector>
#include <string>
#include <iostream>

#include "TFile.h"
#include "TSpline.h"
#include "TLegend.h"

// for local
#include "util.h"
#include "drawset.h"
#include "mathfunc.h"
#include "xrayline.hh"
#include "fitparam.hh"
#include "fitUser.hh"


int fituser(TString target="he4",Int_t nrebin=8)
{
  TString fname="hoge.root";
  TString hname="fuga";
  TString ln="KHe4LAlpha";
  if ( (target=="he3") || (target=="He3") ) {
    fname=Form("%s/hene_phc_sprmon160_301_he3_spline.root",dout.Data());
    hname="hene_phc_kheton160_301_runs";
    ln="KHe3LAlpha";
  }else if ( (target=="he4") || (target=="He4") ) {
    fname=Form("%s/hene_phc_sprmon320_424_he4_spline.root",dout.Data());
    hname="hene_phc_kheton320_424_runs";
    ln="KHe4LAlpha";
  }else{
    fname=Form("%s/hene_phc_sprmon320_424_he4_spline.root",dout.Data());
    hname="hene_phc_kheton320_424_runs";
    ln="FeKa";
  }
  if (gSystem->AccessPathName(fname)) {
    cout << "missing file: " << fname << endl;
    return 1;
  }
  cout << "file name: " << fname << endl;
  cout << "hist name: " << hname << endl;
  cout << "line name: " << ln << endl;
  TFile f(fname,"read");
  TH1F *h1 = (TH1F*)f.Get(hname);
  Double_t pp=0;
  Bool_t LET=true;Bool_t HET=false;Bool_t BG=true;
  Bool_t ENE=true;Bool_t PLOT=true;Bool_t SAVE=true;
  FitUserXray(h1,nrebin,pp,ln,LET,HET,BG,ENE,PLOT,SAVE);
  return 0;
}


#endif
