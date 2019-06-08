#ifndef HISTFUNC_HH
#define HISTFUNC_HH 1

#include "TROOT.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TFile.h"
#include "TGraph.h"
#include "TSpectrum.h"
#include "TSpline.h"
#include "TList.h"
#include "TPolyMarker.h"
#include "TRandom3.h"

#include <vector>
#include <iostream>
#include <fstream>

#include "util.h"
#include "drawset.h"


// peak search for multipule xray peaks
// returns a vector contains peak positions
inline std::vector<double> peaksearch(TH1F *h1, Int_t nrebin=20,Int_t np=6,
                                      Double_t p_th=0.03, Double_t p_sigma=3)
{
  Int_t maxnpeak=TMath::Abs(np);
  TString hname = h1->GetName();
  TString htitle = h1->GetTitle();
  // use a copy to keep the original histogram
  TH1F *h = (TH1F*)h1->Clone(Form("%s_p",hname.Data()));
  h->Rebin(nrebin);
  h->SetTitle(htitle);
  Double_t bw = h->GetBinWidth(1);

  // Use TSpectrum to find the peak candidates
  TSpectrum *s = new TSpectrum(maxnpeak);
  Int_t nfound = s->Search(h,p_sigma,"nodraw",p_th);
  std::vector<double> pp;// peak positions
  Double_t *xpeaks = s->GetPositionX();
  cout << Form("found %d peaks on %s",nfound,hname.Data()) << endl;
  for (Int_t p=0; p<nfound; p++) { pp.push_back(xpeaks[p]); }
  std::sort(pp.begin(), pp.end());
  for (Int_t p=0; p<nfound; p++) {cout << p+1 << ": " << pp[p] << endl;}
  
  // canvas to check peaks by eyes
  SetAxisHist(h,"adc [ch]",Form("Counts / %.2f ch",bw),61,61);
  h->GetXaxis()->SetRangeUser(pp[0]-5000,pp[nfound-1]+5000);
  TCanvas *c_ps = new TCanvas("c_ps","c_ps",10,10,810,610);
  c_ps->Divide(1,2);
  c_ps->cd(2);
  gPad->SetGrid(0,0);
  gPad->SetLogy(1);
  h->DrawCopy("hist");
  c_ps->cd(1);
  gPad->SetGrid(0,0);
  Double_t maxh = h->GetMaximum();
  h->GetYaxis()->SetRangeUser(0.,maxh*1.5);
  h->DrawCopy("hist");
  TLegend *lhist = new TLegend(0.70, 0.40, 0.85, 0.80);
  lhist->SetFillColor(0);
  lhist->SetFillStyle(0);
  lhist->SetTextSize(0.04);
  lhist->AddEntry(h->GetName(),h->GetName(),"pe");
  for (Int_t p=0; p<nfound; p++) {
    lhist->AddEntry((TObject*)0,Form("pp[%d]: %.0f",p,pp[p]),"");  
  }
  lhist->Draw("same");
  TList *functions = h->GetListOfFunctions();
  TPolyMarker *pm = (TPolyMarker*)functions->FindObject("TPolyMarker");
  pm->Draw("same");
  
  // Estimate background using TSpectrum::Background
  TH1 *hb = s->Background(h,20,"same");
  if (hb) c_ps->Update();
  //TString coutname=Form("%s/%s/%s_%s.pdf",fig_dir.Data(),ps_tag.Data(),ps_tag.Data(),hname.Data());
  //c_ps->SaveAs(coutname);
  
  if (hb) hb->Delete();
  pm->Delete();
  functions->Delete();
  c_ps->Close();
  lhist->Delete();
  s->Delete();
  h->Delete();

  return pp;
}


void PeakSearchTest(Int_t nrebin=10,Double_t p_th=0.020,
                    Int_t run=381, TString fname="hoge.root")
{
  // channels you want to check peak search
  std::vector<int> chans{13,15,21,47,53,61,99,169,265,293,427};

  TFile f(fname);
  for (int i=0; i<(int)chans.size(); i++) {
    f.cd();
    TString keyname = Form("hpht_phc_on%d_ch%d",run,chans[i]);
    TH1F *h1 = (TH1F*)f.Get(keyname);
    TString hname = h1->GetName();
    auto pp = peaksearch(h1,nrebin,6,p_th,3);
    const Int_t np = int(pp.size());
    if (np!=6) {
      cout << "number of peaks looks strange " << hname << endl;
      ofstream ofst(badpsfile, ios::app);
      if (!ofst) {cout << "file open error " << badpsfile << endl; exit(1);}
      ofst << getdt() << " " << hname << " " << np << endl;
    }
    h1->Delete();
  }
  f.cd();
  f.Close();
  return;
}

inline TGraph* GetGraph(Int_t np,
                         Double_t *x, Double_t *y,
                         TString hname="hist", TString tag="")
{
  TString g_name = Form("%s_%s",tag.Data(),hname.Data());
  TGraph *g = new TGraph(np,x,y);
  g->SetName(g_name);
  g->SetTitle(g_name);
  return g;
}

inline TGraphErrors* GetGraphErrors(Int_t np,
                                    Double_t *x, Double_t *y,
                                    Double_t *xErr, Double_t *yErr,
                                    TString hname="hist", TString tag="")
{
  TString g_name = Form("%s_%s",tag.Data(),hname.Data());
  TGraphErrors *g = new TGraphErrors(np,x,y,xErr,yErr);
  g->SetName(g_name);
  g->SetTitle(g_name);
  return g;
}


inline TSpline3* GetSpline3(TGraph *g, TString hname="hist", TString tag="")
{
  TString sp3_name = Form("%s_%s",tag.Data(),hname.Data());
  TSpline3 *sp3 = new TSpline3(sp3_name,g);
  sp3->SetName(sp3_name);
  return sp3;
}

inline TGraph* GetSpline3DerivativeGraph(TSpline3 *sp1,
                                         TString hname="hist",
                                         const Int_t npder=1000,
                                         Double_t ini=4000., Double_t fin=10000.)
{
  Double_t x[npder];
  Double_t y[npder];
  Double_t step = (double)(fin - ini)/npder;
  for (int i=0; i<npder; i++) {
    x[i] = ini+step*i;
    y[i] = 1./sp1->Derivative(x[i]);
  }
  TGraph *gsp3der = GetGraph(npder,x,y,hname,sp3meander_tag);
  return gsp3der;
}

inline TGraphErrors* GetGraphMeans(TString fname, TString hname="hist")
{
  TFile f(fname,"read");
  if (f.IsZombie()) {
    cout << "Warning: cannot open file " << fname << endl;
    TGraphErrors *dummy=nullptr;
    return dummy;
  }
  TString g_name = Form("%s_%s",gmean_tag.Data(),hname.Data());
  TGraphErrors *g = (TGraphErrors*)f.Get(g_name);
  f.Close();
  return g;
}

inline TSpline3* GetSpline3(TString fname, TString hname="hist")
{
  TFile f(fname,"read");
  if (f.IsZombie()) {
    cout << "Warning: cannot open file " << fname << endl;
    TSpline3 *dummy=nullptr;
    return dummy;
  }
  TString sp3_name = Form("%s_%s",sp3calib_tag.Data(),hname.Data());
  TSpline3 *sp3 = (TSpline3*)f.Get(sp3_name);
  f.Close();
  return sp3;
}

inline std::vector<double> GetFitMeans(TString fname, TString hname="hist")
{
  TFile f(fname,"read");
  if (f.IsZombie()) {
    cout << "Warning: cannot open file " << fname << endl;
    std::vector<double> dummy;
    return dummy;
  }
  TString gmean_name = Form("%s_%s",gmean_tag.Data(),hname.Data());
  TGraphErrors *g = (TGraphErrors*)f.Get(gmean_name);
  Int_t np = g->GetN();
  Double_t *pp = g->GetY();
  std::vector<double> gm;
  for (int i=0; i<np; i++) {
    gm.push_back(pp[i]);
    cout << Form("graph mean [%d]: %.1f",i,pp[i]) << endl;
  }
  g->Delete();
  f.Close();
  return gm;
}

inline TGraphErrors* GetGraphMeansDir(TString fname, TString hname="hist")
{
  TFile f(fname,"read");
  if (f.IsZombie()) {
    cout << "Warning: cannot open file " << fname << endl;
    TGraphErrors *dummy=nullptr;
    return dummy;
  }
  TDirectory *fd = (TDirectory*)f.GetDirectory(hname);
  TString g_name = Form("%s_%s",gmean_tag.Data(),hname.Data());
  TGraphErrors *g = (TGraphErrors*)fd->Get(g_name);
  f.Close();
  return g;
}

inline TSpline3* GetSpline3MeanDir(TString fname, TString hname="hist")
{
  TFile f(fname,"read");
  if (f.IsZombie()) {
    cout << "Warning: cannot open file " << fname << endl;
    TSpline3 *dummy=nullptr;
    return dummy;
  }
  TDirectory *fd = (TDirectory*)f.GetDirectory(hname);
  TString sp3_name = Form("%s_%s",sp3mean_tag.Data(),hname.Data());
  TSpline3 *sp3 = (TSpline3*)fd->Get(sp3_name);
  f.Close();
  return sp3;
}

inline TSpline3* GetSpline3CalibDir(TString fname, TString hname="hist")
{
  TFile f(fname,"read");
  if (f.IsZombie()) {
    cout << "Warning: cannot open file " << fname << endl;
    TSpline3 *dummy=nullptr;
    return dummy;
  }
  TDirectory *fd = (TDirectory*)f.GetDirectory(hname);
  TString sp3_name = Form("%s_%s",sp3calib_tag.Data(),hname.Data());
  TSpline3 *sp3 = (TSpline3*)fd->Get(sp3_name);
  f.Close();
  return sp3;
}

inline TSpline3* GetSpline3Dir(TString fname, TString hname="hist")
{
  TFile f(fname,"read");
  if (f.IsZombie()) {
    cout << "Warning: cannot open file " << fname << endl;
    TSpline3 *dummy=nullptr;
    return dummy;
  }
  TDirectory *fd = (TDirectory*)f.GetDirectory(hname);
  TString sp3_name = Form("%s_%s",sp3calib_tag.Data(),hname.Data());
  TSpline3 *sp3 = (TSpline3*)fd->Get(sp3_name);
  f.Close();
  return sp3;
}

inline std::vector<double> GetFitMeansDir(TString fname, TString hname="hist")
{
  TFile f(fname,"read");
  if (f.IsZombie()) {
    cout << "Warning: cannot open file " << fname << endl;
    std::vector<double> dummy;
    return dummy;
  }
  TDirectory *fd = (TDirectory*)f.GetDirectory(hname);
  TString gmean_name = Form("%s_%s",gmean_tag.Data(),hname.Data());
  TGraphErrors *g = (TGraphErrors*)fd->Get(gmean_name);
  Int_t np = g->GetN();
  Double_t *pp = g->GetY();
  std::vector<double> gm;
  for (int i=0; i<np; i++) {
    gm.push_back(pp[i]);
    cout << Form("graph mean [%d]: %.1f",i,pp[i]) << endl;
  }
  g->Delete();
  f.Close();
  return gm;
}


inline Double_t GetE2CFromSpline(Double_t ene, TString hname="hist")
{
  TString fname=Form("%s/%s/%s_%s.root",
                     par_dir.Data(),calib_tag.Data(),calib_tag.Data(),hname.Data());
  TFile f(fname,"read");
  if (f.IsZombie()) {
    cout << "Warning: cannot open file " << fname << endl;
    cout << "use default values 0.6 for e2c" << endl;
    return 0.6;
  }
  TString sp3mean_name = Form("%s_%s",sp3mean_tag.Data(),hname.Data());
  TSpline3 *g = (TSpline3*)f.Get(sp3mean_name);
  Double_t e2c = 1./g->Derivative(ene);
  g->Delete();
  f.Close();
  return e2c;
}


inline Double_t GetADCFromSpline(Double_t ene, TString hname="hist")
{
  TString fname=Form("%s/%s/%s_%s.root",
                     par_dir.Data(),calib_tag.Data(),calib_tag.Data(),hname.Data());
  TFile f(fname,"read");
  if (f.IsZombie()) {
    cout << "Warning: cannot open file " << fname << endl;
    return 0.;
  }
  TString sp3mean_name = Form("%s_%s",sp3mean_tag.Data(),hname.Data());
  TSpline3 *g = (TSpline3*)f.Get(sp3mean_name);
  Double_t adc = g->Eval(ene);
  g->Delete();
  f.Close();
  return adc;
}

inline TH1F* GetEnergyConvertedHist(TH1 *h,TSpline3 *sp3calib,
                                    TString hname="hist",TString tag="",
                                    Double_t low=0., Double_t high=20000.)
{
  Int_t enebin=(int)(high-low)*4.;//0.25 eV/ch
  TString henename=Form("%s_%s",tag.Data(),hname.Data());
  delete gROOT->FindObject(henename);
  TH1F *hene = new TH1F(henename,henename,enebin,low,high);
  Int_t nbin=h->GetNbinsX();
  Double_t initch = h->GetXaxis()->GetBinLowEdge(1);
  Double_t lastch = h->GetXaxis()->GetBinLowEdge(nbin);
  for (Int_t j=0; j<(int)(lastch-initch)+1; j++) {
    Double_t chw = h->GetBinContent(j+1);
    Double_t adc = h->GetBinCenter(j+1);
    TRandom3 *rnd = new TRandom3();
    rnd->SetSeed();
    for(Int_t k=0; k<chw; k++) {
      Double_t enesp3 = (sp3calib->Eval(adc+rnd->Rndm()-0.5));
      hene->Fill(enesp3);
    }
    rnd->Delete();
  }
  return hene;
}

inline TH1F* EnergyConversion(TH1F *h1,TString hname_calib,Bool_t ZERO=false)
{
  // this should be run after making calibration curve with TSpline3
  // convert ADC ch into energy for other histogram "h1"
  // hname_calib is the name of "hcalib" used for energy calibration
  TString hname=h1->GetName();
  TString htitle=h1->GetTitle();
  if (contains_check(hname_calib,badch)) return nullptr;
  TString fcalib=Form("%s/%s/%s_%s.root",par_dir.Data(),calib_tag.Data(),calib_tag.Data(),hname_calib.Data());
  if (ZERO) fcalib=Form("%s/%s/%s_%s_zero.root",par_dir.Data(),calib_tag.Data(),calib_tag.Data(),hname_calib.Data());
  if (gSystem->AccessPathName(fcalib)) return nullptr;
  TSpline3 *sp3calib = GetSpline3(fcalib,hname_calib);
  TH1F *hene = GetEnergyConvertedHist(h1,sp3calib,hname,hene_tag);
  sp3calib->Delete();
  return hene;
}

inline TH1F* EnergyConversion(TH1F *h1,TSpline3 *sp3calib)
{
  // this should be run after making calibration curve with TSpline3
  // convert ADC ch into energy for other histogram "h1"
  TString hname=h1->GetName();
  TString htitle=h1->GetTitle();
  if (contains_check(hname,badch)) return nullptr;
  TH1F *hene = GetEnergyConvertedHist(h1,sp3calib,hname,hene_tag);
  return hene;
}





// --------------------------------------------------------
// --------------------------------------------------------








inline TGraphErrors *TFResidue(TF1* func, TGraphErrors* ge,
                               Double_t frll, Double_t frul,Bool_t PLOT=false)
{
  Int_t N = ge->GetN();
  Double_t *x = ge->GetX();
  Double_t *xErr = ge->GetEX();
  Double_t *y = ge->GetY();
  Double_t *yErr = ge->GetEY();
  Double_t fy[N],res[N],resErr[N];
  
  for (Int_t i=0; i<N; i++) {
    fy[i] = func->Eval(x[i]);
    res[i] = y[i]-fy[i];
    resErr[i] = yErr[i];
  }
  
  TGraphErrors *gres = new TGraphErrors(N,x,res,xErr,resErr);
  gres->SetTitle("fit residue");
  
  if (PLOT) {
    Double_t upper=2.0;
    Double_t lower=-2.0;
    Double_t upperch = func->Eval(frul);
    Double_t lowerch = func->Eval(frll);
    
    gROOT->SetStyle("Plain");
    TCanvas *cresi = new TCanvas("cresi","residue graph",1,1,800,600);
    gPad->SetGrid();//gPad->SetLogy();
    gPad->SetTicky(1);
    //gStyle->SetOptFit(1111);
    //gStyle->SetFitFormat("14.7g");
    gStyle->SetOptStat(kFALSE);
    SetAxisMarkerTGraphErrors(gres,"voltage [V]","residue val - f(x) [ch]",62,62,20,1);
    gres->GetXaxis()->SetRangeUser(frll,frul);
    gres->GetYaxis()->SetRangeUser(lower,upper);
    gres->Draw("AP");

    TF1 *f0 = new TF1("f0","x",lowerch,upperch); 
    f0->SetNpx(10000);
    TGaxis *A1 = new TGaxis(frll,upper,frul,upper,"f0",510,"-"); 
    A1->SetTitle("ADC channel"); 
    A1->Draw();

    cresi->Update();
    
    //if(SAVE) {
    //cresi->Print(Form("pdf/"+type+"_residue_sdd%d.pdf",chid));
    //cout << Form("pdf/"+type+"_residue_sdd%d.pdf",chid) << " was created." << endl;
    //}
  }
  
  return gres;
}





inline TGraphErrors *TFResidueY(TF1* func, TGraphErrors* ge,
                                Double_t frll, Double_t frul,Bool_t PLOT=false)
{
  Int_t N = ge->GetN();
  Double_t *x = ge->GetX();
  Double_t *xErr = ge->GetEX();
  Double_t *y = ge->GetY();
  Double_t *yErr = ge->GetEY();
  Double_t fy[N],res[N],resErr[N];
  
  for (Int_t i=0; i<N; i++) {
    fy[i] = func->GetX(y[i]);
    res[i] = x[i]-fy[i];
    resErr[i] = xErr[i];// <- now zero !
  }
  
  TGraphErrors *gres = new TGraphErrors(N,y,res,yErr,resErr);
  gres->SetTitle("fit residue");
  
  if (PLOT) {
    Double_t upper=2.0;
    Double_t lower=-2.0;
    Double_t upperch = func->Eval(frul);
    Double_t lowerch = func->Eval(frll);
    
    gROOT->SetStyle("Plain");
    TCanvas *cresi = new TCanvas("cresi","residue graph",1,1,800,600);
    gPad->SetGrid();//gPad->SetLogy();
    gPad->SetTicky(1);
    //gStyle->SetOptFit(1111);
    //gStyle->SetFitFormat("14.7g");
    gStyle->SetOptStat(kFALSE);
    
    SetAxisMarkerTGraphErrors(gres,"voltage [V]","residue val - f(x) [ch]",62,62,20,1);
    gres->GetXaxis()->SetRangeUser(frll,frul);
    gres->GetYaxis()->SetRangeUser(lower,upper);
    gres->Draw("AP");

    TF1 *f0 = new TF1("f0","x",lowerch,upperch); 
    f0->SetNpx(10000);
    TGaxis *A1 = new TGaxis(frll,upper,frul,upper,"f0",510,"-"); 
    A1->SetTitle("ADC channel"); 
    A1->Draw();

    cresi->Update();
    
    //if(SAVE) {
    //cresi->Print(Form("pdf/"+type+"_residue_sdd%d.pdf",chid));
    //cout << Form("pdf/"+type+"_residue_sdd%d.pdf",chid) << " was created." << endl;
    //}
  }
  
  return gres;
}



inline TGraphErrors *GetHistResidue(TF1* func, TH1F* h,
                                    Double_t frll, Double_t frul,
                                    Bool_t PLOT=false)
{
  Int_t lowBin = h->FindBin(frll);
  Int_t upBin = h->FindBin(frul);
  Int_t nBin = upBin-lowBin+1;

  Double_t x[nBin];
  Double_t y[nBin];
  Double_t fy[nBin];
  Double_t res[nBin];
  Double_t resError[nBin];

  for (Int_t i=0; i<nBin; i++) {
    x[i] = h->GetBinCenter(i+lowBin);
    y[i] = h->GetBinContent(i+lowBin);
    fy[i] = func->Eval(x[i]);
    res[i] = y[i]-fy[i];
    resError[i] = h->GetBinError(i+lowBin);
  }

  TGraphErrors *gres = new TGraphErrors(nBin,x,res,0,resError);
  gres->SetTitle("fit residue");
  
  if (PLOT) {
    gROOT->SetStyle("Plain");
    TCanvas *cresi = new TCanvas("cresi","residue graph",1,1,800,600);
    gPad->SetGrid();
    gPad->SetTicks();
    gStyle->SetOptStat(kFALSE);
    SetAxisMarkerTGraphErrors(gres,"Energy (eV)","residue val(x) - f(x)",62,62,24,1);
    gres->Draw("AP");
    cresi->Update();
  }

  return gres;
}



inline TGraph *GetHistResiduePoint(TF1* func, TH1F* h,
                                   Double_t frll, Double_t frul,
                                   Bool_t PLOT=false)
{
  Int_t lowBin = h->FindBin(frll);
  Int_t upBin = h->FindBin(frul);
  Int_t nBin = upBin-lowBin+1;

  Double_t x[nBin];
  Double_t y[nBin];
  Double_t fy[nBin];
  Double_t res[nBin];
  Double_t resError[nBin];

  for (Int_t i=0; i<nBin; i++) {
    x[i] = h->GetBinCenter(i+lowBin);
    y[i] = h->GetBinContent(i+lowBin);
    fy[i] = func->Eval(x[i]);
    res[i] = y[i]-fy[i];
    resError[i] = h->GetBinError(i+lowBin);
  }

  TGraph *gres = new TGraph(nBin,x,res);
  gres->SetTitle("fit residue");
  
  if (PLOT) {
    gROOT->SetStyle("Plain");
    TCanvas *cresi = new TCanvas("cresi","residue graph",1,1,800,600);
    gPad->SetGrid();
    gPad->SetTicks();
    gStyle->SetOptStat(kFALSE);
    SetAxisMarkerTGraph(gres,"Energy (eV)","residue val(x) - f(x)",62,62,24,1);
    gres->Draw("AP");
    cresi->Update();
  }

  return gres;
}



//TGraph *GetHistResidueTwoSigmaLinePlus(TF1* func, TH1F* h,
inline TGraph *GetHistResidueTwoSigmaLinePlus(TH1F* h,Double_t frll, Double_t frul,
                                              Bool_t PLOT=false)
{
  Int_t lowBin = h->FindBin(frll);
  Int_t upBin = h->FindBin(frul);
  Int_t nBin = upBin-lowBin+1;
  
  Double_t x[nBin];
  Double_t resError[nBin];

  for (Int_t i=0; i<nBin; i++) {
    x[i] = h->GetBinCenter(i+lowBin);
    resError[i] = +2.0*(h->GetBinError(i+lowBin));
  }

  TGraph *g = new TGraph(nBin,x,resError);
  g->SetTitle("fit residue");
  
  if (PLOT) {
    gROOT->SetStyle("Plain");
    TCanvas *c15 = new TCanvas("c15","residue 2sigma plue line",1,1,800,600);
    gPad->SetGrid();
    gPad->SetTicks();
    gStyle->SetOptStat(kFALSE);
    SetAxisMarkerTGraph(g,"Energy (eV)","residue val(x) - f(x)",62,62,24,1);
    g->Draw("APC");
    c15->Update();
  }

  return g;
}



//TGraph *GetHistResidueTwoSigmaLineMinus(TF1* func, TH1F* h,
inline TGraph *GetHistResidueTwoSigmaLineMinus(TH1F* h,Double_t frll, Double_t frul,
                                                Bool_t PLOT=false)
{
  Int_t lowBin = h->FindBin(frll);
  Int_t upBin = h->FindBin(frul);
  Int_t nBin = upBin-lowBin+1;
  
  Double_t x[nBin];
  Double_t resError[nBin];

  for (Int_t i=0; i<nBin; i++) {
    x[i] = h->GetBinCenter(i+lowBin);
    resError[i] = -2.0*(h->GetBinError(i+lowBin));
  }

  TGraph *g = new TGraph(nBin,x,resError);
  g->SetTitle("fit residue");
  
  if (PLOT) {
    gROOT->SetStyle("Plain");
    TCanvas *c25 = new TCanvas("c25","residue 2sigma minus line",1,1,800,600);
    gPad->SetGrid();
    gPad->SetTicks();
    gStyle->SetOptStat(kFALSE);
    SetAxisMarkerTGraph(g,"Energy (eV)","residue val(x) - f(x)",62,62,24,1);
    g->Draw("APC");
    c25->Update();
  }

  return g;
}



//TGraph *GetHistResidueTwoSigmaLinePlus(TF1* func, TH1F* h,
inline TGraph *GetHistResidueThreeSigmaLinePlus(TH1F* h,Double_t frll, Double_t frul,
                                                Bool_t PLOT=false)
{
  Int_t lowBin = h->FindBin(frll);
  Int_t upBin = h->FindBin(frul);
  Int_t nBin = upBin-lowBin+1;
  
  Double_t x[nBin];
  Double_t resError[nBin];

  for (Int_t i=0; i<nBin; i++) {
    x[i] = h->GetBinCenter(i+lowBin);
    resError[i] = +3.0*(h->GetBinError(i+lowBin));
  }

  TGraph *g = new TGraph(nBin,x,resError);
  g->SetTitle("fit residue");
  
  if (PLOT) {
    gROOT->SetStyle("Plain");
    TCanvas *c15 = new TCanvas("c15","residue 2sigma plue line",1,1,800,600);
    gPad->SetGrid();
    gPad->SetTicks();
    gStyle->SetOptStat(kFALSE);
    SetAxisMarkerTGraph(g,"Energy (eV)","residue val(x) - f(x)",62,62,24,1);
    g->Draw("APC");
    c15->Update();
  }

  return g;
}



//TGraph *GetHistResidueTwoSigmaLineMinus(TF1* func, TH1F* h,
inline TGraph *GetHistResidueThreeSigmaLineMinus(TH1F* h,Double_t frll, Double_t frul,
                                                 Bool_t PLOT=false)
{
  Int_t lowBin = h->FindBin(frll);
  Int_t upBin = h->FindBin(frul);
  Int_t nBin = upBin-lowBin+1;
  
  Double_t x[nBin];
  Double_t resError[nBin];

  for (Int_t i=0; i<nBin; i++) {
    x[i] = h->GetBinCenter(i+lowBin);
    resError[i] = -3.0*(h->GetBinError(i+lowBin));
  }

  TGraph *g = new TGraph(nBin,x,resError);
  g->SetTitle("fit residue");
  
  if (PLOT) {
    gROOT->SetStyle("Plain");
    TCanvas *c25 = new TCanvas("c25","residue 2sigma minus line",1,1,800,600);
    gPad->SetGrid();
    gPad->SetTicks();
    gStyle->SetOptStat(kFALSE);
    SetAxisMarkerTGraph(g,"Energy (eV)","residue val(x) - f(x)",62,62,24,1);
    g->Draw("APC");
    c25->Update();
  }

  return g;
}



inline TCanvas* DrawResidueGraphs(TH1F *h0,TF1 *func,Double_t low,Double_t high,Bool_t PLOT=false)
{
  // residuals
  TGraph *res = GetHistResiduePoint(func,h0,low,high,PLOT);
  res->SetName("residue");
  TGraph *gp = GetHistResidueTwoSigmaLinePlus(h0,low,high,PLOT);
  gp->SetName("plus2sigma");
  TGraph *gm = GetHistResidueTwoSigmaLineMinus(h0,low,high,PLOT);

  TCanvas *cres = new TCanvas("cres","residue",1,1,800,600);
  gPad->SetGrid();
  gPad->SetTicks();
  gStyle->SetOptStat(kFALSE);
  res->SetTitle("fit residue");
  SetAxisMarkerTGraph(res,"Energy (eV)","residue val(x) - f(x)",62,62,7,1);
  res->Draw("APL");
  gp->SetLineWidth(2);
  gm->SetLineWidth(2);
  gp->Draw("PL");
  gm->Draw("PL");
  cres->Update();

  return cres;
}




inline TGraphErrors *GetGraphResidue(TF1* func, TGraphErrors* ge,
                                     Bool_t PLOT=false)
{
  Int_t npoint   = ge->GetN();
  Double_t *x    = ge->GetX();
  Double_t *y    = ge->GetY();
  //Double_t *xErr = ge->GetEX();// unused
  Double_t *yErr = ge->GetEY();

  Double_t fy[npoint];
  Double_t res[npoint];
  Double_t resError[npoint];

  for (Int_t i=0; i<npoint; i++) {
    fy[i]  = func->Eval(x[i]);
    res[i] = y[i]-fy[i];
    resError[i] = yErr[i];
    cout << Form("x:%f, y:%f, fy:%f, res:%f",x[i],y[i],fy[i],res[i]) << endl;
  }

  TGraphErrors *gres = new TGraphErrors(npoint,x,res,0,resError);
  gres->SetTitle("fit residue");
  
  if (PLOT) {
    gROOT->SetStyle("Plain");
    TCanvas *cresi = new TCanvas("cresi","residue graph",1,1,800,600);
    gPad->SetGrid();
    gPad->SetTicks();
    gStyle->SetOptStat(kFALSE);
    SetAxisMarkerTGraphErrors(gres,"Energy (eV)","residue val(x) - f(x)",62,62,24,1);
    gres->Draw("APL");
    cresi->Update();
  }

  return gres;
}


#endif
