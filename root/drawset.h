#ifndef _DRAWSET_H
#define _DRAWSET_H 1

#include "TMath.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TMinuit.h"
#include "TCanvas.h"

// about font number: see pp.132 of ROOT Guide
// "Graphics and the Graphical User Interface"
// "62" is recomended

inline Bool_t SetColorWidthDraw(TF1 *func, Int_t color, Int_t width, TString opt)
{
  func->SetNpx(10000); func->SetLineColor(color); 
  func->SetLineWidth(width); func->Draw(opt);
  return true;
}

inline Bool_t SetAxisHist(TH1F *h, TString Xtitle, TString Ytitle, Int_t Xfont, Int_t Yfont)
{
  h->GetXaxis()->SetTitle(Xtitle); 
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->SetTitle(Ytitle);  
  h->GetYaxis()->CenterTitle();
  h->GetXaxis()->SetLabelFont(Xfont);
  h->GetYaxis()->SetLabelFont(Yfont);
  h->GetXaxis()->SetTitleFont(Xfont);
  h->GetYaxis()->SetTitleFont(Yfont);
  h->GetXaxis()->SetTitleOffset(1.1);
  h->GetYaxis()->SetTitleOffset(1.3);
  return true;
}

inline Bool_t SetAxisHist2D(TH2F *h, TString Xtitle, TString Ytitle, Int_t Xfont, Int_t Yfont)
{
  h->GetXaxis()->SetTitle(Xtitle); 
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->SetTitle(Ytitle);  
  h->GetYaxis()->CenterTitle();
  h->GetXaxis()->SetLabelFont(Xfont);
  h->GetYaxis()->SetLabelFont(Yfont);
  h->GetXaxis()->SetTitleFont(Xfont);
  h->GetYaxis()->SetTitleFont(Yfont);
  h->GetXaxis()->SetTitleOffset(1.1);
  h->GetYaxis()->SetTitleOffset(1.3);
  return true;
}

inline Bool_t SetAxisTF1(TF1 *f, TString Xtitle, TString Ytitle, Int_t Xfont, Int_t Yfont)
{
  f->GetXaxis()->SetTitle(Xtitle); 
  f->GetXaxis()->CenterTitle();
  f->GetYaxis()->SetTitle(Ytitle);  
  f->GetYaxis()->CenterTitle();
  f->GetXaxis()->SetLabelFont(Xfont);
  f->GetYaxis()->SetLabelFont(Yfont);
  f->GetXaxis()->SetTitleFont(Xfont);
  f->GetYaxis()->SetTitleFont(Yfont);
  return true;
}

inline void SetAxisTGraph(TGraph *f, TString Xtitle, TString Ytitle, Int_t Xfont, Int_t Yfont)
{
  f->GetXaxis()->SetTitle(Xtitle);
  f->GetXaxis()->CenterTitle();
  f->GetYaxis()->SetTitle(Ytitle);
  f->GetYaxis()->CenterTitle();
  f->GetXaxis()->SetLabelFont(Xfont);
  f->GetYaxis()->SetLabelFont(Yfont);
  f->GetXaxis()->SetTitleFont(Xfont);
  f->GetYaxis()->SetTitleFont(Yfont);
  return;
}

inline Bool_t SetAxisTGraphErrors(TGraphErrors *f, TString Xtitle, TString Ytitle, Int_t Xfont, Int_t Yfont)
{
  f->GetXaxis()->SetTitle(Xtitle); 
  f->GetXaxis()->CenterTitle();
  f->GetYaxis()->SetTitle(Ytitle);  
  f->GetYaxis()->CenterTitle();
  f->GetXaxis()->SetLabelFont(Xfont);
  f->GetYaxis()->SetLabelFont(Yfont);
  f->GetXaxis()->SetTitleFont(Xfont);
  f->GetYaxis()->SetTitleFont(Yfont);
  return true;
}

inline Bool_t SetAxisTMultiGraph(TMultiGraph *f, TString Xtitle, TString Ytitle, Int_t Xfont, Int_t Yfont)
{
  f->GetXaxis()->SetTitle(Xtitle); 
  f->GetXaxis()->CenterTitle();
  f->GetYaxis()->SetTitle(Ytitle);  
  f->GetYaxis()->CenterTitle();
  f->GetXaxis()->SetLabelFont(Xfont);
  f->GetYaxis()->SetLabelFont(Yfont);
  f->GetXaxis()->SetTitleFont(Xfont);
  f->GetYaxis()->SetTitleFont(Yfont);
  return true;
}

inline Bool_t SetAxisMarkerTGraph(TGraph *f, TString Xtitle, TString Ytitle, Int_t Xfont, Int_t Yfont, Int_t mstyle, Int_t mcolor)
{
  f->GetXaxis()->SetTitle(Xtitle); 
  f->GetXaxis()->CenterTitle();
  f->GetYaxis()->SetTitle(Ytitle);  
  f->GetYaxis()->CenterTitle();
  f->GetXaxis()->SetLabelFont(Xfont);
  f->GetYaxis()->SetLabelFont(Yfont);
  f->GetXaxis()->SetTitleFont(Xfont);
  f->GetYaxis()->SetTitleFont(Yfont);
  f->SetMarkerStyle(mstyle);
  f->SetMarkerColor(mcolor);
  f->SetLineColor(mcolor);
  return true;
}

inline Bool_t SetMarkerTGraphErrors(TGraphErrors *f, TString title="",Int_t mstyle=20,Int_t mlcolor=1, Int_t fstyle=0)
{

  f->SetTitle(title);
  f->SetMarkerStyle(mstyle);
  f->SetMarkerColor(mlcolor);
  f->SetLineColor(mlcolor);
  f->SetFillStyle(fstyle);
  return true;
}

inline Bool_t SetAxisMarkerTGraphErrors(TGraphErrors *f, TString Xtitle, TString Ytitle, Int_t Xfont, Int_t Yfont, Int_t mstyle, Int_t mcolor)
{
  f->GetXaxis()->SetTitle(Xtitle); 
  f->GetXaxis()->CenterTitle();
  f->GetYaxis()->SetTitle(Ytitle);  
  f->GetYaxis()->CenterTitle();
  f->GetXaxis()->SetLabelFont(Xfont);
  f->GetYaxis()->SetLabelFont(Yfont);
  f->GetXaxis()->SetTitleFont(Xfont);
  f->GetYaxis()->SetTitleFont(Yfont);
  f->SetMarkerStyle(mstyle);
  f->SetMarkerColor(mcolor);
  f->SetLineColor(mcolor);
  return true;
}

inline Bool_t ShowChisquareNDF(TF1 *func)
{
  Double_t chisq = func->GetChisquare();
  Int_t    ndf   = func->GetNDF();
  cout << Form("chisq/ndf = %3.3f",chisq/ndf) << endl;
  return true;
}


/*
inline TCanvas *DrawMinuCont(Int_t npoints, Int_t par1, Int_t par2,
			     TString Xtitle="", TString Ytitle="", 
			     Int_t Xfont=132, Int_t Yfont=132,
			     Bool_t PLOT = true, Bool_t SAVE=false, 
			     TString fname = "hoge.pdf")
{
  TCanvas *cont = new TCanvas("cont","contours",750,700);
  cont->cd();
  gMinuit->SetErrorDef(9);
  TGraph *gr3 = (TGraph*)gMinuit->Contour(npoints,par1,par2);
  gr3->SetLineColor(1);
  gr3->SetFillColor(0);
  
  gMinuit->SetErrorDef(4);
  TGraph *gr2 = (TGraph*)gMinuit->Contour(npoints,par1,par2);
  gr2->SetLineColor(1);
  gr2->SetFillColor(0);
  
  gMinuit->SetErrorDef(1);
  TGraph *gr1 = (TGraph*)gMinuit->Contour(npoints,par1,par2);
  gr1->SetLineColor(1);
  gr1->SetFillColor(0);
  gr3->SetTitle("fit contour (1,2 and 3 #sigma)");
  gr3->GetYaxis()->SetLabelSize(0.03);
  gr3->GetYaxis()->SetTitleSize(0.03);
  gr3->GetYaxis()->SetTitleOffset(1.65);
  gr3->GetXaxis()->SetLabelSize(0.03);
  gr3->GetXaxis()->SetTitleSize(0.03);

  SetAxisTGraph(gr3,Xtitle,Ytitle,Xfont,Yfont);
  SetAxisTGraph(gr2,Xtitle,Ytitle,Xfont,Yfont);
  SetAxisTGraph(gr1,Xtitle,Ytitle,Xfont,Yfont);

  if (PLOT==true) {
    gr3->Draw("ac");
    gr2->Draw("c");
    gr1->Draw("c");
  }
  cont->Update();

  if (SAVE==true) {
    cont->Print("./fig/"+fname);
    cont->Print("./tmp_cont.C");
    //cout << "./figure/cont/"+fname+" was created" << endl;
  }
  return cont;
}
*/

#endif
