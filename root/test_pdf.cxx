#include "TRandom.h"
#include "TCanvas.h"
#include "TLegend.h"
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
  TH1F *htailsm = new TH1F("htailsm","htailsm",nbin,xlow,xhigh);
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
  TH1F *htailsm = new TH1F("htailsm","htailsm",nbin,xlow,xhigh);
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
