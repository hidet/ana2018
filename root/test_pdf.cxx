#include <iostream>
#include <fstream>

#include "TRandom.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TF1Convolution.h"

#include "drawset.h"
#include "xrayline.hh"



void test_pdf(TString line_name="CoKa", Double_t fwhm=6.0)
{
  xrayline XL = get_xl(line_name);
  Double_t ampl=XL.normalized_intensities[0];
  Double_t mean=XL.energies[0];
  Double_t sigma=fwhm/2.355;
  
  Double_t xlow=mean*0.95;
  Double_t xhigh=mean*1.02;
  TString pdfname=Form("%s.pdf(x)",line_name.Data());
  TF1 *lrnz = new TF1("lrnz",pdfname,xlow,xhigh);
  lrnz->SetNpx(10000);
  
  const int npar=6;
  Double_t tail_ratio=0.2;
  Double_t tail_beta=50.;// eV
  TF1 *voigt_tail = new TF1("voigt_tail",XL,xlow,xhigh,npar);
  Double_t par[npar]={1.0,ampl,mean,sigma,tail_ratio,tail_beta};
  voigt_tail->SetParameters(par);
  TString paranames[npar]={"eV per ch","Voigt ampl","Voigt mean","Voigt sigma","LE tail ratio","LE tail slope"};
  for (Int_t i=0; i<npar; i++) {voigt_tail->SetParName(i,paranames[i].Data());}
  voigt_tail->SetNpx(10000);
  
  Int_t nbin = int(xhigh-xlow);
  TH1F *hsim = new TH1F("hsim","hsim",nbin,xlow,xhigh);
  TH1F *hvsim = new TH1F("hvsim","hvsim",nbin,xlow,xhigh);
  TH1F *hl = new TH1F("hl","hl",nbin,xlow,xhigh);
  TH1F *htail = new TH1F("htail","htail",nbin,xlow,xhigh);
  //TH1F *htailsm = new TH1F("htailsm","htailsm",nbin,xlow,xhigh);
  Double_t bw = hl->GetBinWidth(1);
      
  for (int i=0;i<1e6;i++) {
    Double_t l = lrnz->GetRandom();
    hl->Fill(l);
    l += gRandom->Gaus(0.,sigma);
  }

  // Lorenztian -> Exp tail
  TString ftail=Form("exp(x/%f)*((x < 0) ? 1 : 0)",tail_beta);
  TF1 *fexp = new TF1("fexp",ftail,-tail_beta*10.,0);
  for (int i=1;i<nbin+1;i++) {
    Double_t c = hl->GetBinContent(i);
    Double_t m = hl->GetBinCenter(i);
    for (int j=0;j<c;j++) {
      Double_t t = fexp->GetRandom();
      htail->Fill(m+t);
    }
  }
  fexp->Delete();

  // add Lorenztian and exp tail with considering the ratio
  hsim->Add(hl,htail,1.-tail_ratio,tail_ratio);

  // Gauss smearing
  for (int i=1;i<nbin+1;i++) {
    Double_t c = hsim->GetBinContent(i);
    Double_t m = hsim->GetBinCenter(i);
    for (int j=0;j<c;j++) {
      Double_t t = gRandom->Gaus(0.,sigma);
      hvsim->Fill(m+t);
    }
  }

  voigt_tail->FixParameter(0, 1.0);
  voigt_tail->SetParLimits(1, 0., 1e7);
  voigt_tail->SetParLimits(2, mean-5., mean+5.);
  voigt_tail->SetParLimits(3, 1., 10.);
  voigt_tail->SetParLimits(4, 0.05, 0.3);
  voigt_tail->SetParLimits(5, 1., 60.);
  hvsim->Fit(voigt_tail,"ELR0");

  
  TCanvas *c1 = new TCanvas("c1","c1",1);
  gPad->SetGrid(0,0);
  gPad->SetTicks();
  gPad->SetLogy();
  SetAxisHist(hl,"Energy [eV]",Form("Counts / %.2f eV",bw),61,61);
  hl->SetLineColor(4);
  hl->SetTitle(Form("LE Tail smearing test, fwhm=%.1f, LE ratio=%.2f, LE slope=%.1f",fwhm,tail_ratio,tail_beta));
  hl->Draw("hist");

  //lrnz->Draw("same");
  
  htail->SetLineColor(6);
  htail->Draw("histsame");
  
  hvsim->Draw("histsame");
  voigt_tail->SetLineColor(3);
  voigt_tail->Draw("same");

  TLegend *lhist = new TLegend(0.2, 0.6, 0.6, 0.85);
  lhist->SetFillColor(0);
  lhist->SetFillStyle(0);
  lhist->SetTextSize(0.03);
  lhist->AddEntry(hl->GetName(),"Lorentzian","e");
  lhist->AddEntry(htail->GetName(),"Exp tail","e");
  lhist->AddEntry(hvsim->GetName(),"Simulated Voigt + LE Tail","e");
  lhist->AddEntry(voigt_tail->GetName(),"Tatsuno's fit func","l");
  lhist->Draw();

  c1->Update();
  c1->SaveAs(Form("./fig/test_pdf_%s.pdf",line_name.Data()));
  
}









// try to use TF1convolution
void draw_pdf_conv(TString line_name="CoKa", Double_t fwhm=6.0)
{
  xrayline XL = get_xl(line_name);
  Double_t ampl=XL.normalized_intensities[0];
  Double_t mean=XL.energies[0];
  Double_t sigma=fwhm/2.355;
  
  Double_t xlow=mean*0.95;
  Double_t xhigh=mean*1.02;
  TString pdfname=Form("%s.pdf(x)",line_name.Data());
  TF1 *lrnz = new TF1("lrnz",pdfname,xlow,xhigh);
  lrnz->SetNpx(10000);
  
  const int npar=6;
  Double_t tail_ratio=0.2;
  Double_t tail_beta=50.;// eV
  TF1 *voigt_tail = new TF1("voigt_tail",XL,xlow,xhigh,npar);
  Double_t par[npar]={1.0,ampl,mean,sigma,tail_ratio,tail_beta};
  voigt_tail->SetParameters(par);
  TString paranames[npar]={"eV per ch","Voigt ampl","Voigt mean","Voigt sigma","LE tail ratio","LE tail slope"};
  for (Int_t i=0; i<npar; i++) {voigt_tail->SetParName(i,paranames[i].Data());}
  voigt_tail->SetNpx(10000);
  
  Int_t nbin = int(xhigh-xlow);
  TH1F *hsim = new TH1F("hsim","hsim",nbin,xlow,xhigh);
  TH1F *hvsim = new TH1F("hvsim","hvsim",nbin,xlow,xhigh);
  TH1F *hl = new TH1F("hl","hl",nbin,xlow,xhigh);
  TH1F *htail = new TH1F("htail","htail",nbin,xlow,xhigh);
  //TH1F *htailsm = new TH1F("htailsm","htailsm",nbin,xlow,xhigh);
  Double_t bw = hl->GetBinWidth(1);
      
  for (int i=0;i<1e6;i++) {
    Double_t l = lrnz->GetRandom();
    hl->Fill(l);
    l += gRandom->Gaus(0.,sigma);
  }

  TString ftail=Form("exp(x/%f)*((x < 0) ? 1 : 0)",tail_beta);
  TF1 *fexp = new TF1("fexp",ftail,-tail_beta*10.,0);
  TF1Convolution *conv = new TF1Convolution(lrnz,fexp);
  conv->SetRange(xlow*0.8,xhigh*1.2);
  conv->SetNofPointsFFT(10000);
  TF1 *fconv = new TF1("fconv",*conv,xlow,xhigh,conv->GetNpar());
  fconv->SetNpx(10000);
  for (int i=0;i<1e6;i++) {
    Double_t t = fconv->GetRandom();
    htail->Fill(t);
  }
 
  // add Lorenztian and exp tail with considering the ratio
  hsim->Add(hl,htail,1.-tail_ratio,tail_ratio);

  // Gauss smearing
  for (int i=1;i<nbin+1;i++) {
    Double_t c = hsim->GetBinContent(i);
    Double_t m = hsim->GetBinCenter(i);
    for (int j=0;j<c;j++) {
      Double_t t = gRandom->Gaus(0.,sigma);
      hvsim->Fill(m+t);
    }
  }
  
  voigt_tail->FixParameter(0, 1.0);
  voigt_tail->SetParLimits(1, 0., 1e7);
  voigt_tail->SetParLimits(2, mean-5., mean+5.);
  voigt_tail->SetParLimits(3, 1., 10.);
  voigt_tail->SetParLimits(4, 0.05, 0.3);
  voigt_tail->SetParLimits(5, 1., 60.);
  hvsim->Fit(voigt_tail,"ELR0");
  
  TCanvas *c1 = new TCanvas("c1","c1",1);
  gPad->SetGrid(0,0);
  gPad->SetTicks();
  gPad->SetLogy();
  SetAxisHist(hl,"Energy [eV]",Form("Counts / %.2f eV",bw),61,61);
  hl->SetLineColor(4);
  hl->SetTitle(Form("LE Tail smearing test, fwhm=%.1f, LE ratio=%.2f, LE slope=%.1f",fwhm,tail_ratio,tail_beta));
  hl->Draw("hist");
  //hl->DrawNormalized("hist");

  //lrnz->Draw("same");
  
  htail->SetLineColor(6);
  htail->Draw("histsame");
  //htail->DrawNormalized("histsame");
  
  hvsim->Draw("histsame");
  voigt_tail->SetLineColor(3);
  voigt_tail->Draw("same");

  TLegend *lhist = new TLegend(0.2, 0.6, 0.6, 0.85);
  lhist->SetFillColor(0);
  lhist->SetFillStyle(0);
  lhist->SetTextSize(0.03);
  lhist->AddEntry(hl->GetName(),"Lorentzian","e");
  lhist->AddEntry(htail->GetName(),"Exp tail","e");
  lhist->AddEntry(hvsim->GetName(),"Simulated Voigt + LE Tail","e");
  lhist->AddEntry(voigt_tail->GetName(),"New fit func","l");
  lhist->Draw();

  c1->Update();
  //c1->SaveAs(Form("./fig/test_pdf_%s.pdf",line_name.Data()));

}




TH1F* get_hist(TString line_name="TiKa",Double_t fwhm=6.0,Int_t nev=1e6,
               Double_t xlow=0., Double_t xhigh=10000., Int_t nbin=10000,
               Bool_t LETAIL=true, Bool_t SAVE=false, TString foutname="hoge.dat")
{
  xrayline XL = get_xl(line_name);
  Double_t ampl=XL.normalized_intensities[0];
  Double_t mean=XL.energies[0];
  Double_t sigma=fwhm/2.355;
  
  //Double_t xlow=mean*0.95;
  //Double_t xhigh=mean*1.02;
  TString pdfname=Form("%s.pdf(x)",line_name.Data());
  TF1 *lrnz = new TF1("lrnz",pdfname,xlow,xhigh);
  lrnz->SetNpx(10000);
  
  const int npar=6;
  Double_t tail_ratio=0.2;
  if (!LETAIL) tail_ratio=0.0;
  Double_t tail_beta=50.;// eV
  TF1 *voigt_tail = new TF1("voigt_tail",XL,xlow,xhigh,npar);
  Double_t par[npar]={1.0,ampl,mean,sigma,tail_ratio,tail_beta};
  voigt_tail->SetParameters(par);
  TString paranames[npar]={"eV per ch","Voigt ampl","Voigt mean","Voigt sigma","LE tail ratio","LE tail slope"};
  for (Int_t i=0; i<npar; i++) {voigt_tail->SetParName(i,paranames[i].Data());}
  voigt_tail->SetNpx(10000);
  
  TH1F *hsim = new TH1F("hsim","hsim",nbin,xlow,xhigh);
  TH1F *hvsim = new TH1F(Form("hvsim_%s",line_name.Data()),Form("hvsim_%s",line_name.Data()),nbin,xlow,xhigh);
  TH1F *hl = new TH1F("hl","hl",nbin,xlow,xhigh);
  TH1F *htail = new TH1F("htail","htail",nbin,xlow,xhigh);
  //TH1F *htailsm = new TH1F("htailsm","htailsm",nbin,xlow,xhigh);
  Double_t bw = hl->GetBinWidth(1);
      
  for (int i=0;i<nev;i++) {
    Double_t l = lrnz->GetRandom();
    hl->Fill(l);
    l += gRandom->Gaus(0.,sigma);
  }

  // Lorenztian -> Exp tail
  TString ftail=Form("exp(x/%f)*((x < 0) ? 1 : 0)",tail_beta);
  TF1 *fexp = new TF1("fexp",ftail,-tail_beta*10.,0);
  for (int i=1;i<nbin+1;i++) {
    Double_t c = hl->GetBinContent(i);
    Double_t m = hl->GetBinCenter(i);
    for (int j=0;j<c;j++) {
      Double_t t = fexp->GetRandom();
      htail->Fill(m+t);
    }
  }
  fexp->Delete();

  // add Lorenztian and exp tail with considering the ratio
  hsim->Add(hl,htail,1.-tail_ratio,tail_ratio);

  ofstream ofst;
  if (SAVE) {
    ofst.open(foutname, ios::app);
    if (!ofst) {cout << "file open error " << endl; exit(1);}
  }
  // Gauss smearing
  for (int i=1;i<nbin+1;i++) {
    Double_t c = hsim->GetBinContent(i);
    Double_t m = hsim->GetBinCenter(i);
    for (int j=0;j<c;j++) {
      Double_t t = gRandom->Gaus(0.,sigma);
      hvsim->Fill(m+t);
      if (SAVE) ofst << m+t << endl;
    }
  }
   
  hsim->Delete();
  hl->Delete();
  htail->Delete();
  //voigt_tail->FixParameter(0, 1.0);
  //voigt_tail->SetParLimits(1, 0., 1e7);
  //voigt_tail->SetParLimits(2, mean-5., mean+5.);
  //voigt_tail->SetParLimits(3, 1., 10.);
  //voigt_tail->SetParLimits(4, 0.05, 0.3);
  //voigt_tail->SetParLimits(5, 1., 60.);
  //hvsim->Fit(voigt_tail,"ELR0");

  return hvsim;
}




void TestCalib(Int_t nka=1e4, Bool_t LETAIL=true)
{
  Int_t nkb=0.15*nka;
  Double_t xlow=0.;
  Double_t xhigh=10000.;
  Int_t nbin=10000;

  TString frootname="test_calib_withtail.root";
  TString foutname="test_calib_withtail.dat";
  
  if (!LETAIL) {
    frootname="test_calib.root";
    foutname="test_calib.dat";
  }
  TFile f(frootname,"recreate");

  TH1F *hcrka = get_hist("CrKa",5.6,nka,xlow,xhigh,nbin,LETAIL,true,foutname);
  TH1F *hcrkb = get_hist("CrKb",5.6,nkb,xlow,xhigh,nbin,LETAIL,true,foutname);
  TH1F *hcoka = get_hist("CoKa",6.0,nka,xlow,xhigh,nbin,LETAIL,true,foutname);
  TH1F *hcokb = get_hist("CoKb",7.0,nkb,xlow,xhigh,nbin,LETAIL,true,foutname);
  TH1F *hcuka = get_hist("CuKa",8.0,nka,xlow,xhigh,nbin,LETAIL,true,foutname);

  TList *list = new TList;
  list->Add(hcrka);
  list->Add(hcrkb);
  list->Add(hcoka);
  list->Add(hcokb);
  list->Add(hcuka);

  TH1F *hsum = (TH1F*)hcrka->Clone("hsum");
  hsum->Reset();
  hsum->Merge(list);

  hsum->Write();
  f.Close();

  return;
}



void hoge(Int_t nka=1e4)
{
  Int_t nkb=0.15*nka;
  Double_t xlow=0.;
  Double_t xhigh=10000.;
  Int_t nbin=10000;
  
  TH1F *htika = get_hist("TiKa",5.3,nka,xlow,xhigh,nbin);
  TH1F *htikb = get_hist("TiKb",5.3,nkb,xlow,xhigh,nbin);
  TH1F *hvka  = get_hist("VKa",5.4,nka,xlow,xhigh,nbin);
  TH1F *hvkb  = get_hist("VKb",5.4,nkb,xlow,xhigh,nbin);
  TH1F *hcrka = get_hist("CrKa",5.6,nka,xlow,xhigh,nbin);
  TH1F *hcrkb = get_hist("CrKb",5.6,nkb,xlow,xhigh,nbin);
  TH1F *hmnka = get_hist("MnKa",6.0,nka,xlow,xhigh,nbin);
  TH1F *hmnkb = get_hist("MnKb",6.0,nkb,xlow,xhigh,nbin);
  TH1F *hfeka = get_hist("FeKa",6.5,nka,xlow,xhigh,nbin);
  TH1F *hfekb = get_hist("FeKb",6.5,nkb,xlow,xhigh,nbin);
  TH1F *hcoka = get_hist("CoKa",7.0,nka,xlow,xhigh,nbin);
  TH1F *hcokb = get_hist("CoKb",7.0,nkb,xlow,xhigh,nbin);
  TH1F *hnika = get_hist("NiKa",7.5,nka,xlow,xhigh,nbin);
  TH1F *hnikb = get_hist("NiKb",7.5,nkb,xlow,xhigh,nbin);
  TH1F *hcuka = get_hist("CuKa",8.0,nka,xlow,xhigh,nbin);
  TH1F *hcukb = get_hist("CuKb",8.0,nkb,xlow,xhigh,nbin);
  

  TList *list = new TList;
  list->Add(htika);
  list->Add(htikb);
  list->Add(hvka);
  list->Add(hvkb);
  list->Add(hcrka);
  list->Add(hcrkb);
  list->Add(hmnka);
  list->Add(hmnkb);
  list->Add(hfeka);
  list->Add(hfekb);
  list->Add(hcoka);
  list->Add(hcokb);
  list->Add(hnika);
  list->Add(hnikb);
  list->Add(hcuka);
  list->Add(hcukb);
  
  TH1F *hsum = (TH1F*)htika->Clone("hsum");
  hsum->Reset();
  hsum->Merge(list);

  //hsum->Draw("hist");

  htika -> SetLineColor(2);
  htikb -> SetLineColor(2);
  hvka  -> SetLineColor(3);
  hvkb  -> SetLineColor(3);
  hcrka -> SetLineColor(4);
  hcrkb -> SetLineColor(4);
  hmnka -> SetLineColor(5);
  hmnkb -> SetLineColor(5);
  hfeka -> SetLineColor(6);
  hfekb -> SetLineColor(6);
  hcoka -> SetLineColor(7);
  hcokb -> SetLineColor(7);
  hnika -> SetLineColor(8);
  hnikb -> SetLineColor(8);
  hcuka -> SetLineColor(9);
  hcukb -> SetLineColor(9);
  
 // htika -> Draw("histsame");
 // htikb -> Draw("histsame");
  hvka  -> Draw("histsame");
  hvkb  -> Draw("histsame");
  //hcrka -> Draw("histsame");
  //hcrkb -> Draw("histsame");
  hmnka -> Draw("histsame");
  hmnkb -> Draw("histsame");
  //hfeka -> Draw("histsame");
  //hfekb -> Draw("histsame");
  hcoka -> Draw("histsame");
  hcokb -> Draw("histsame");
  //hnika -> Draw("histsame");
  //hnikb -> Draw("histsame");
  hcuka -> Draw("histsame");
  hcukb -> Draw("histsame");
  
  
}
