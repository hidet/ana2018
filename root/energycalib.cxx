#ifndef ENERGYCALIB_CXX
#define ENERGYCALIB_CXX 1

#include <vector>
#include <iostream>

#include "TSystem.h"
#include "TFile.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TRandom3.h"

#include "util.h"
#include "drawset.h"
#include "histfunc.hh"
#include "fitElem.hh"

void menu()
{
  cout << "-- Functions --" << endl;
  cout << "FitMultiLines()" <<endl;
  cout << "FitMultiLinesEne()" <<endl;
  cout << "LoopChs()" <<endl;
  cout << "energycalib()" <<endl;
  cout << endl;
  cout << "-- How to use --" << endl;
  cout << "run -q energycalib.cxx" << endl;
  // run -q 'energycalib.cxx("sprmon","he4")'
  return;
}

inline void FitMultiLines(TH1F *h1,Int_t rebin=4,
                          Bool_t LET=true,Bool_t HET=false,Bool_t BG=false,
                          Bool_t ENE=false,Bool_t PLOT=false,Bool_t SAVE=true,
                          TFile *fout=nullptr)
{
  // h1: histogram to be fitted
  // fout: TFile to save results
  TString hname=h1->GetName();
  TString htitle=h1->GetTitle();
  Bool_t ITR=false;// iterative for KBeta parameters
  std::vector<double> pp;// peak positions
  TSpline3 *sp3=nullptr;
  TString fcalib=Form("%s/%s/%s_%s.root",
                      par_dir.Data(),calib_tag.Data(),calib_tag.Data(),hname.Data());
  if (!gSystem->AccessPathName(fcalib)) {
    ITR=true;
    TFile f(fcalib,"read");
    TString gmean_name = Form("%s_%s",gmean_tag.Data(),hname.Data());
    TGraphErrors *g = (TGraphErrors*)f.Get(gmean_name);
    Int_t np = g->GetN();
    Double_t *p = g->GetY();
    for (int i=0; i<np; i++) {
      pp.push_back(p[i]);
      cout << Form("fit mean [%d]: %.1f",i,p[i]) << endl;
    }
    g->Delete();
    sp3 = (TSpline3*)f.Get(Form("%s_%s",sp3mean_tag.Data(),hname.Data()));
    f.Close();
  } else {
    pp = peaksearch(h1,10,6,0.020,3);
  }
  const Int_t np = int(pp.size());
  Double_t ene[np],  eneErr[np];
  Double_t mean[np], meanErr[np];
  const Int_t np_zero = int(pp.size())+1;
  Double_t ene_zero[np_zero],  eneErr_zero[np_zero];
  Double_t mean_zero[np_zero], meanErr_zero[np_zero];
  ene_zero[0]=0.;   eneErr_zero[0]=0.;
  mean_zero[0]=0.;  meanErr_zero[0]=0.;
  std::vector<fitparam> fps{};
  if (np==5) {
    TString peaknames[5]={"CrKa","CrKb","CoKa","CoKb","CuKa"};
    for(Int_t i=0; i != int(pp.size()); ++i) {
      if (ITR) fps.push_back(FitElemXrayItr(h1,sp3,rebin,peaknames[i].Data(),LET,HET,BG,ENE,PLOT,SAVE));
      else fps.push_back(FitElemXray(h1,rebin,pp[i],peaknames[i].Data(),LET,HET,BG,ENE,PLOT,SAVE));
      mean[i]    = fps[i].mean;
      meanErr[i] = fps[i].meanErr;
      ene[i]     = fps[i].ene;
      eneErr[i]  = fps[i].eneErr;
      mean_zero[i+1]    = fps[i].mean;
      meanErr_zero[i+1] = fps[i].meanErr;
      ene_zero[i+1]     = fps[i].ene;
      eneErr_zero[i+1]  = fps[i].eneErr;
    }
  } else if (np==6) {
    TString peaknames[6]={"CrKa","CrKb","CoKa","CoKb","CuKa","CuKb"};
    for(Int_t i=0; i != int(pp.size()); ++i) {
      if (ITR) fps.push_back(FitElemXrayItr(h1,sp3,rebin,peaknames[i].Data(),LET,HET,BG,ENE,PLOT,SAVE));
      else fps.push_back(FitElemXray(h1,rebin,pp[i],peaknames[i].Data(),LET,HET,BG,ENE,PLOT,SAVE));
      mean[i] = fps[i].mean;
      meanErr[i] = fps[i].meanErr;
      ene[i] =  fps[i].ene;
      eneErr[i] =  fps[i].eneErr;
      mean_zero[i+1]    = fps[i].mean;
      meanErr_zero[i+1] = fps[i].meanErr;
      ene_zero[i+1]     = fps[i].ene;
      eneErr_zero[i+1]  = fps[i].eneErr;
    }
  }else{
    cout << "number of peaks looks strange " << hname << endl;
    ofstream ofst(badpsfile, ios::app);
    if (!ofst) {cout << "file open error " << badpsfile << endl; exit(1);}
    ofst << getdt() << " " << hname << " " << np << endl;
    if (sp3) sp3->Delete();
    return;
  }
  if (sp3) sp3->Delete();
  if (ENE) return;
  //------------------------------------------
  // energy calibration with TSpline3
  TGraphErrors *gmean = GetGraphErrors(np,ene,mean,eneErr,meanErr,hname,gmean_tag);
  TSpline3 *sp3mean = GetSpline3(gmean,hname,sp3mean_tag);
  TGraph *gsp3meander = GetSpline3DerivativeGraph(sp3mean,hname,1000,4000.,10000.);
  TGraphErrors *gcalib = GetGraphErrors(np,mean,ene,meanErr,eneErr,hname,gcalib_tag);
  TSpline3 *sp3calib = GetSpline3(gcalib,hname,sp3calib_tag);
  TH1F *hene = GetEnergyConvertedHist(h1,sp3calib,hname,hene_tag);
  //------------------------------------------
  TString hname_zero=hname+"_zero";
  TGraphErrors *gmean_zero  = GetGraphErrors(np_zero,ene_zero,mean_zero,eneErr_zero,meanErr_zero,hname_zero,gmean_tag);
  TSpline3 *sp3mean_zero    = GetSpline3(gmean_zero,hname_zero,sp3mean_tag);
  TGraph *gsp3meander_zero  = GetSpline3DerivativeGraph(sp3mean_zero,hname_zero,2000,0.,10000.);
  TGraphErrors *gcalib_zero = GetGraphErrors(np_zero,mean_zero,ene_zero,meanErr_zero,eneErr_zero,hname_zero,gcalib_tag);
  TSpline3 *sp3calib_zero   = GetSpline3(gcalib_zero,hname_zero,sp3calib_tag);
  TH1F *hene_zero           = GetEnergyConvertedHist(h1,sp3calib_zero,hname_zero,hene_tag);
  //------------------------------------------
  TString calibname=Form("%s_%s",calib_tag.Data(),hname.Data());
  TCanvas *ccalib = new TCanvas(calibname,calibname,1,1,600,800);
  ccalib->Divide(1,2);
  ccalib->cd(1);
  gPad->SetGrid(0,0);
  gPad->SetTicks();
  SetAxisTGraphErrors(gmean,"Energy (eV)","adc (ch)",61,61);
  gmean->GetXaxis()->SetLimits(ene[0]-1000,ene[np-1]+1000);
  gmean->SetMarkerStyle(20);
  gmean->Draw("AP");
  sp3mean->SetLineColor(4);
  sp3mean->Draw("same");
  TLegend *lmean = new TLegend(0.15, 0.65, 0.5, 0.80);
  lmean->SetFillColor(0);
  lmean->SetFillStyle(0);
  lmean->SetTextSize(0.03);
  lmean->AddEntry(gmean->GetName(),"fitted means","ep");
  lmean->AddEntry(sp3mean->GetName(),"cubic spline interpolate","l");
  lmean->Draw("same");
  ccalib->cd(2);
  gPad->SetGrid(0,0);
  gPad->SetTicks();
  SetAxisTGraph(gsp3meander,"Energy (eV)","eV/ch",61,61);
  gsp3meander->GetXaxis()->SetLimits(ene[0]-1000,ene[np-1]+1000);
  gsp3meander->SetMinimum(0.45);
  gsp3meander->SetMaximum(0.90);
  gsp3meander->SetLineColor(4);
  gsp3meander->Draw("APC");
  TLegend *le2c = new TLegend(0.15, 0.65, 0.5, 0.80);
  le2c->SetFillColor(0);
  le2c->SetFillStyle(0);
  le2c->SetTextSize(0.04);
  le2c->AddEntry(gsp3meander->GetName(),"e2c from 1./cubic_spline_derivative","l");
  le2c->Draw("same");
  ccalib->Update();
  //------------------------------------------
  TString calibname_zero=Form("%s_%s",calib_tag.Data(),hname_zero.Data());
  TCanvas *ccalib_zero = new TCanvas(calibname_zero,calibname_zero,1,1,600,800);
  ccalib_zero->Divide(1,2);
  ccalib_zero->cd(1);
  gPad->SetGrid(0,0);
  gPad->SetTicks();
  SetAxisTGraphErrors(gmean_zero,"Energy (eV)","adc (ch)",61,61);
  gmean_zero->GetXaxis()->SetLimits(ene_zero[0]-1000,ene_zero[np_zero-1]+1000);
  gmean_zero->SetMarkerStyle(20);
  gmean_zero->Draw("AP");
  sp3mean_zero->SetLineColor(4);
  sp3mean_zero->Draw("same");
  TLegend *lmean_zero = new TLegend(0.15, 0.65, 0.5, 0.80);
  lmean_zero->SetFillColor(0);
  lmean_zero->SetFillStyle(0);
  lmean_zero->SetTextSize(0.03);
  lmean_zero->AddEntry(gmean_zero->GetName(),"fitted means","ep");
  lmean_zero->AddEntry(sp3mean_zero->GetName(),"cubic spline interpolate","l");
  lmean_zero->Draw("same");
  ccalib_zero->cd(2);
  gPad->SetGrid(0,0);
  gPad->SetTicks();
  SetAxisTGraph(gsp3meander_zero,"Energy (eV)","eV/ch",61,61);
  gsp3meander_zero->GetXaxis()->SetLimits(ene_zero[0]-1000,ene_zero[np_zero-1]+1000);
  gsp3meander_zero->SetMinimum(0.45);
  gsp3meander_zero->SetMaximum(0.90);
  gsp3meander_zero->SetLineColor(4);
  gsp3meander_zero->Draw("APC");
  TLegend *le2c_zero = new TLegend(0.15, 0.65, 0.5, 0.80);
  le2c_zero->SetFillColor(0);
  le2c_zero->SetFillStyle(0);
  le2c_zero->SetTextSize(0.04);
  le2c_zero->AddEntry(gsp3meander_zero->GetName(),"e2c from 1./cubic_spline_derivative","l");
  le2c_zero->Draw("same");
  ccalib_zero->Update();
  //------------------------------------------
  if (SAVE) {
    TString coutmean=Form("%s/%s/%s.pdf",fig_dir.Data(),calib_tag.Data(),calibname.Data());
    ccalib->SaveAs(coutmean);
    TString coutmean_zero=Form("%s/%s/%s.pdf",fig_dir.Data(),calib_tag.Data(),calibname_zero.Data());
    ccalib_zero->SaveAs(coutmean_zero);
    TH1F *h = (TH1F*)h1->Clone(Form("%s_c",hname.Data()));
    h->SetTitle(htitle);
    //------------------------------------------
    if (fout) {// for merge channels
      fout->cd();
      h           ->Write();
      hene        ->Write();
      gmean       ->Write();
      sp3mean     ->Write();
      gsp3meander ->Write();
      gcalib      ->Write();
      sp3calib    ->Write();
      ccalib      ->Write();
      hene_zero        ->Write();
      gmean_zero       ->Write();
      sp3mean_zero     ->Write();
      gsp3meander_zero ->Write();
      gcalib_zero      ->Write();
      sp3calib_zero    ->Write();
      ccalib_zero      ->Write();
    }
    TString foutname=Form("%s/%s/%s_%s.root",
                          par_dir.Data(),calib_tag.Data(),calib_tag.Data(),hname.Data());
    TFile foutch(foutname,"recreate");
    foutch.cd();
    h           ->Write();
    hene        ->Write();
    gmean       ->Write();
    sp3mean     ->Write();
    gsp3meander ->Write();
    gcalib      ->Write();
    sp3calib    ->Write();
    ccalib      ->Write();
    TString foutname_zero=Form("%s/%s/%s_%s.root",
                               par_dir.Data(),calib_tag.Data(),calib_tag.Data(),hname_zero.Data());
    TFile foutch_zero(foutname_zero,"recreate");
    foutch_zero.cd();
    h           ->Write();
    hene_zero        ->Write();
    gmean_zero       ->Write();
    sp3mean_zero     ->Write();
    gsp3meander_zero ->Write();
    gcalib_zero      ->Write();
    sp3calib_zero    ->Write();
    ccalib_zero      ->Write();
    //------------------------------------------
    h           ->Delete();
    hene        ->Delete();
    hene_zero   ->Delete();
    gmean       ->Delete();
    sp3mean     ->Delete();
    gsp3meander ->Delete();
    gcalib      ->Delete();
    sp3calib    ->Delete();
    ccalib->Close();
    gmean_zero       ->Delete();
    sp3mean_zero     ->Delete();
    gsp3meander_zero ->Delete();
    gcalib_zero      ->Delete();
    sp3calib_zero    ->Delete();
    ccalib_zero->Close();
    foutch.Close();
    foutch_zero.Close();
    //------------------------------------------
    cout << foutname << " has been created" << endl;
    cout << foutname_zero << " has been created" << endl;
  }
  return;
}

inline void FitMultiLinesEne(TH1F *h1,Int_t rebin=4,
                             Bool_t LET=true,Bool_t HET=false,Bool_t BG=false,
                             Bool_t ENE=true,Bool_t PLOT=true,Bool_t SAVE=true)
{
  TString hname=h1->GetName();
  TString htitle=h1->GetTitle();
  if (contains_check(hname,badch)) return;
  Bool_t ITR=false;// iterative for KBeta parameters
  std::vector<double> pp;// peak positions
  TSpline3 *sp3=nullptr;
  TString fcalib=Form("%s/%s/%s_%s.root",
                      par_dir.Data(),calib_tag.Data(),calib_tag.Data(),hname.Data());
  if (!gSystem->AccessPathName(fcalib)) {
    ITR=true;
    TFile f(fcalib,"read");
    TString gmean_name = Form("%s_%s",gmean_tag.Data(),hname.Data());
    TGraphErrors *g = (TGraphErrors*)f.Get(gmean_name);
    Int_t np = g->GetN();
    Double_t *p = g->GetY();
    for (int i=0; i<np; i++) {
      pp.push_back(p[i]);
      cout << Form("fit mean [%d]: %.1f",i,p[i]) << endl;
    }
    g->Delete();
    sp3 = (TSpline3*)f.Get(Form("%s_%s",sp3mean_tag.Data(),hname.Data()));
    f.Close();
  } else {
    pp = peaksearch(h1,10,6,0.020,3);
  }
  std::vector<fitparam> fps{};
  const Int_t np = int(pp.size());
  if (np==5) {
    TString peaknames[5]={"CrKa","CrKb","CoKa","CoKb","CuKa"};
    for(Int_t i=0; i != int(pp.size()); ++i) {
      if (ITR) fps.push_back(FitElemXrayItr(h1,sp3,rebin,peaknames[i].Data(),LET,HET,BG,ENE,PLOT,SAVE));
      else fps.push_back(FitElemXray(h1,rebin,pp[i],peaknames[i].Data(),LET,HET,BG,ENE,PLOT,SAVE));
    }
  } else if (np==6) {
    TString peaknames[6]={"CrKa","CrKb","CoKa","CoKb","CuKa","CuKb"};
    for(Int_t i=0; i != int(pp.size()); ++i) {
      if (ITR) fps.push_back(FitElemXrayItr(h1,sp3,rebin,peaknames[i].Data(),LET,HET,BG,ENE,PLOT,SAVE));
      else fps.push_back(FitElemXray(h1,rebin,pp[i],peaknames[i].Data(),LET,HET,BG,ENE,PLOT,SAVE));
    }
  }else{
    cout << "number of peaks looks strange " << hname << endl;
    ofstream ofst(badpsfile, ios::app);
    if (!ofst) {cout << "file open error " << badpsfile << endl; exit(1);}
    ofst << getdt() << " " << hname << " " << np << endl;
    if (sp3) sp3->Delete();
    return;
  }
  if (sp3) sp3->Delete();
  return;
}

inline void LoopChs(TString hntag="sprmon",Int_t run=381,
                    TString fname="hoge.root",
                    Int_t rebin=4,
                    Bool_t LET=true,Bool_t HET=false,Bool_t BG=false,
                    Bool_t ENE=false,Bool_t PLOT=false,Bool_t SAVE=true)
{
  if (gSystem->AccessPathName(fname)) {
    cout << "Error: cannot open file " << fname << endl;
    exit(1);
  }
  TString foutname=Form("%s/%s/%s/%s_%s%d.root",
                        par_dir.Data(),calib_tag.Data(),hntag.Data(),
                        calib_tag.Data(),hntag.Data(),run);
  backup_rootfile(foutname);
  TFile fout(foutname,"recreate");
  TFile f(fname,"read");
  TList *list = f.GetListOfKeys();
  TObjLink *lnk = list->FirstLink();
  TString checkname=Form("%s_%s%d_ch",hpht_phc.Data(),hntag.Data(),run);
  while (lnk) {
    TString keyname = lnk->GetObject()->GetName();
    if (keyname.Contains(checkname)) {
      f.cd();
      TH1F *h1 = (TH1F*)f.Get(keyname);
      FitMultiLines(h1,rebin,LET,HET,BG,ENE,PLOT,SAVE,&fout);
      h1->Delete();
    }
    lnk=lnk->Next();
  }
  f.cd();
  f.Close();
  fout.cd();
  fout.Close();
  return;
}

int energycalib(TString hntag="sprmon",TString target="hoge")
{
  gROOT->SetBatch(true);
  TString fname="";
  Int_t nruns=1;
  std::vector<int> runs;
  if (target=="he4") {
    fname = fhe4;
    nruns = nruns_he4;
    for (int i=0; i<nruns; ++i) runs.push_back(runs_he4[i]);
  } else if (target=="he3") {
    fname = fhe3;
    nruns = nruns_he3;
    for (int i=0; i<nruns; ++i) runs.push_back(runs_he3[i]);
  } else {// user define e.g.,
    fname = fhe4;
    runs={395};
    //runs={385,386,387,389};
    nruns=int(runs.size());
  }
  if (gSystem->AccessPathName(fname)) {
    cout << "Error: cannot open file " << fname << endl;
    return 1;
  }
  Int_t rebin=4;
  Bool_t LET=true;Bool_t HET=false;Bool_t BG=false;
  Bool_t ENE=false;Bool_t PLOT=false;Bool_t SAVE=true;
  cout << "fname: " << fname << endl;
  cout << "runs: ";
  for (int i=0; i<nruns; ++i) cout << runs[i] << " ";
  cout << endl;
  check_continue();
  for (int i=0;i<nruns;++i) {
    LoopChs(hntag,runs[i],fname,rebin,LET,HET,BG,ENE,PLOT,SAVE);
  }
  return 0;
}

#endif
