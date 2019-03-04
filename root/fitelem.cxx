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
#include "fitElem.hh"


int fitelem(TString ln="FeKa",TString target="he4",TString hntag="kheton",Int_t nrebin=8)
{
  TString fname=Form("%s/hene_phc_sprmon320_424_he4_spline.root",dout.Data());
  TString hname="hene_phc_kheton320_424_runs";
  if ( (target=="he3") || (target=="He3") ) {
    fname=Form("%s/hene_phc_sprmon160_301_he3_spline.root",dout.Data());
    hname=Form("hene_phc_%s160_301_runs",hntag.Data());
  }else if ( (target=="he4") || (target=="He4") ) {
    fname=Form("%s/hene_phc_sprmon320_424_he4_spline.root",dout.Data());
    hname=Form("hene_phc_%s320_424_runs",hntag.Data());
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
  FitElemXray(h1,nrebin,pp,ln,LET,HET,BG,ENE,PLOT,SAVE);
  return 0;
}


#endif
