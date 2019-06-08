#ifndef ENECALIB_HH
#define ENECALIB_HH 1

#include <vector>
#include <string>
#include <iostream>

#include "TH1.h"
#include "TF1.h"
#include "TSystem.h"
#include "TFile.h"
#include "TSpline.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TGraphErrors.h"

// for local
#include "util.h"
#include "drawset.h"
#include "mathfunc.h"
#include "histfunc.hh"
#include "xrayline.hh"
#include "fitparam.hh"
#include "fitelem.hh"

class EneCalib
{
public:
  EneCalib(){}
  EneCalib(TH1F *h1, Int_t nrebin,
           Bool_t LET=true, Bool_t HET=false, Bool_t BG=false) {
    _hname = h1->GetName();
    _hname_zero = _hname+"_zero";
    _hraw = (TH1F*)h1->Clone(_hname+"_c");// clone
    _hraw->SetTitle(h1->GetTitle());
    _nrebin = nrebin;
    _LET=LET;
    _HET=HET;
    _BG=BG;

    // fit parameters
    _paranames.clear();
    _paranames.push_back("eV per ch");
    _paranames.push_back("Voigt ampl");
    _paranames.push_back("Voigt mean");
    _paranames.push_back("Voigt sigma");
    if (_LET) {
      _paranames.push_back("LE tail ratio");
      _paranames.push_back("LE tail beta");
    }
    if (_HET) {
      _paranames.push_back("HE tail ratio");
      _paranames.push_back("HE tail beta");
    }
    if (_BG) {
      _paranames.push_back("BG const");
      _paranames.push_back("BG slope");
    }
  }
  ~EneCalib() {
    if (_fout) _fout->Close();
  }

  void SetPreviousCalib(TString fprename="hoge.root") {
    if (gSystem->AccessPathName(fprename)) {
      cout << "There is no previous calibration file " << fprename << endl;
      return;
    }
    cout << "previous calib: " << fprename << endl;
    _gmean    = GetGraphMeansDir(fprename,_hname);
    _sp3mean  = GetSpline3MeanDir(fprename,_hname);
    _sp3calib = GetSpline3CalibDir(fprename,_hname);
    return;
  }
  
  void OpenOutputFile(TString foutname="hoge.root", TString opt="update") {
    _foutname=foutname;
    Ssiz_t ind=_hname.Index("_ch");
    _tdrname =  _hname(0,ind);
    cout << "hist name: " << _hname << endl;
    //cout << _tdrname << endl;
    if (_fout==nullptr) {
      _fout = new TFile(_foutname,opt);
      //_fout->mkdir(_tdrname);
      //_fout->cd(_tdrname);
      _fout->mkdir(_hname);
      _fout->cd(_hname);
    }
  }

  void CloseOutputFile() {
    if (_fout) _fout->Close();
  }
  
  void SetFitFunctions(vector<TString> lines) {
    _ffs.clear();
    for (Int_t i=0; i<int(lines.size()); i++) {
      // for each xrayline
      xrayline xl = get_xl(lines[i]);
      TString element = xl.GetElement();
      TString linetype = xl.GetLineType();
      TString funcname = Form("%s_%s",fitElem_tag.Data(),_hname.Data());
      TString xl_tag = Form("%s%s",element.Data(),linetype.Data());
      // range
      Double_t mean = xl.energies[0];
      Double_t low  = mean*0.990;
      Double_t high = mean*1.009;
      // fit func
      Int_t npx=1000;
      FitElem ff(xl,_LET,_HET,_BG,funcname,xl_tag,low,high,npx);
      _ffs.push_back(ff);
    }
  }
  void DoPeakSearch(Int_t nrebin, Int_t np, Double_t p_th, Double_t p_sigma,Bool_t SAVE);
  fitparam DoFitRaw(FitElem ff, Double_t pp, TString ln, Bool_t PLOT, Bool_t SAVE, Bool_t QUIET);
  void DoCalibE62(Bool_t PSSAVE,Bool_t FITPLOT,Bool_t FITSAVE,Bool_t QUIET,Bool_t CALSAVE);

  TH1F* GetHistRaw(){return _hraw;}
  TH1F* GetHistEne(){return _hene;}
  TH1F* GetHistEneZero(){return _hene_zero;}

  void Clear() {
    _gmean            =nullptr;
    _sp3mean          =nullptr;
    _gsp3meander      =nullptr;
    _gcalib           =nullptr;
    _sp3calib         =nullptr;
    _gmean_zero       =nullptr;
    _sp3mean_zero     =nullptr;
    _gsp3meander_zero =nullptr;
    _gcalib_zero      =nullptr;
    _sp3calib_zero    =nullptr;
  }
  
  
private:
  TFile *_fout=nullptr;
  TH1F *_hraw,*_hene,*_hene_zero;
  TString _hname,_hname_zero,_htitle,_foutname,_tdrname;
  Int_t _nrebin;
  Bool_t _LET,_HET,_BG,_PLOT;  
  std::vector<TString> _lines{};// calibration lines
  std::vector<FitElem> _ffs{};// fit functions
  std::vector<double> _pp{};// peak positions
  std::vector<TString> _paranames{};// fit parameters
  
  TGraphErrors *_gmean=nullptr,  *_gmean_zero=nullptr;  
  TSpline3 *_sp3mean=nullptr,    *_sp3mean_zero=nullptr;    
  TGraph *_gsp3meander=nullptr,  *_gsp3meander_zero=nullptr;  
  TGraphErrors *_gcalib=nullptr, *_gcalib_zero=nullptr; 
  TSpline3 *_sp3calib=nullptr,   *_sp3calib_zero=nullptr;  
};



void EneCalib::DoPeakSearch(Int_t nrebin=10, Int_t np=6, Double_t p_th=0.020,
                            Double_t p_sigma=3, Bool_t SAVE=false)
{
  // peak search for multipule xray peaks
  // returns a vector contains peak positions
  Int_t maxnpeak=TMath::Abs(np);
  // use a copy to keep the original histogram
  TH1F *h = (TH1F*)_hraw->Clone(Form("%s_p",_hname.Data()));
  h->Rebin(nrebin);
  h->SetTitle(_htitle+"_peack_search");
  Double_t bw = h->GetBinWidth(1);
  // Use TSpectrum to find the peak candidates
  TSpectrum *s = new TSpectrum(maxnpeak);
  Int_t nfound = s->Search(h,p_sigma,"nodraw",p_th);
  Double_t *xpeaks = s->GetPositionX();
  cout << Form("found %d peaks on %s",nfound,_hname.Data()) << endl;
  _pp.clear();
  for (Int_t p=0; p<nfound; p++) { _pp.push_back(xpeaks[p]); }
  std::sort(_pp.begin(), _pp.end());
  for (Int_t p=0; p<nfound; p++) {cout << p+1 << ": " << _pp[p] << endl;}

  if (SAVE) {
    if (_fout) _fout->cd(_hname);
    // canvas to check peaks by eyes
    SetAxisHist(h,"adc [ch]",Form("Counts / %.2f ch",bw),61,61);
    h->GetXaxis()->SetRangeUser(_pp[0]-5000,_pp[nfound-1]+5000);
    TString cname=Form("c_ps_%s",_hname.Data());
    TCanvas *c_ps = new TCanvas(cname,cname,10,10,810,610);
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
      lhist->AddEntry((TObject*)0,Form("pp[%d]: %.0f",p,_pp[p]),"");  
    }
    lhist->Draw("same");
    TList *functions = h->GetListOfFunctions();
    TPolyMarker *pm = (TPolyMarker*)functions->FindObject("TPolyMarker");
    pm->Draw("same");
    // Estimate background using TSpectrum::Background
    TH1 *hb = (TH1*)s->Background(h,20,"same");
    if (hb) c_ps->Update();
    //TString coutname=Form("%s/%s/%s_%s.pdf",fig_dir.Data(),ps_tag.Data(),ps_tag.Data(),_hname.Data());
    //TString coutname=Form("./output/hoge_peacksearch_%s.pdf",_hname.Data());
    //c_ps->SaveAs(coutname);
    c_ps->Write();

    if (hb) hb->Delete();
    pm->Delete();
    functions->Delete();
    lhist->Delete();
    c_ps->Close();
  }
    
  s->Delete();
  h->Delete();
  return;
}


fitparam EneCalib::DoFitRaw(FitElem ff, Double_t pp=0, TString ln="CrKAlpha",
                            Bool_t PLOT=false, Bool_t SAVE=false, Bool_t QUIET=false) {
  // nrebin: rebin to fit for the histogram
  // pp: peak position
  // ln: name of the xray line (e.g., MnKAlpha, MnKa, mnka)
    
  // use a copy to keep the original histogram                                                               
  TH1F *h = (TH1F*)_hraw->Clone(Form("%s_f",_hname.Data()));
  h->Rebin(_nrebin);
  h->SetTitle(_htitle);
  Double_t bw=h->GetBinWidth(1);
    
  // parameters
  Double_t ampl=1000.;
  Double_t e2c=0.6;
  Double_t mean=pp;
  // xrayline
  xrayline xl = ff.GetXrayLine();
  //xl.Print();
  TString element = xl.GetElement();
  TString linetype = xl.GetLineType();
  if (_sp3mean) {
    e2c = 1./_sp3mean->Derivative(xl.energies[0]);
    mean = _sp3mean->Eval(xl.energies[0]);
  }
  Double_t gaussian_sigma=5;
  // if you need to fix the tail parameters...
  Double_t fixed_tail_ratio=0.20;
  Double_t fixed_tail_slope=30.;
  // initial parameters
  std::vector<double> tp{e2c,ampl,mean,gaussian_sigma};
  Int_t ndp=Int_t(tp.size());
  std::vector<double> tpErr(ndp,0.);// 4 elements with 0.
  if (_LET) {
    tp.push_back(0.1);// LE tail ratio
    tp.push_back(50.);// LE tail beta
    tpErr.push_back(0.);// LE tail ratio
    tpErr.push_back(0.);// LE tail beta
  }
  if (_HET) {
    tp.push_back(0.1);// HE tail ratio
    tp.push_back(50.);// HE tail beta
    tpErr.push_back(0.);// HE tail ratio
    tpErr.push_back(0.);// HE tail beta
  }
  if (_BG) {
    tp.push_back(0.);// BG const
    tp.push_back(0.);// BG slope
    tpErr.push_back(0.);// BG const
    tpErr.push_back(0.);// BG slope
  }
  // range
  h->SetAxisRange(mean*0.98,mean*1.02,"X");
  Double_t low = mean*0.99;
  Double_t high = mean*1.01;
  if ((ln=="FeKa") || (ln=="FeKAlpha")) {
    low=mean*0.986;
    high=mean*1.003;
  }
  // fit func
  Int_t npar = ff.GetNpar();
  TString funcname=ff.GetFuncName();
  TF1 func(funcname,ff,low,high,npar);
  func.SetNpx(ff.GetNpX());
  func.SetParameters(&tp.front());
  func.SetParErrors(&tpErr.front());
  for (Int_t i=0; i<npar; i++) {func.SetParName(i,_paranames[i].Data());}
  // fit parameter boundary
  func.SetParLimits(0, e2c*0.5, e2c*1.5);
  if (linetype=="KBeta") func.FixParameter(0, e2c);
  func.SetParLimits(1, 10., 1e8);
  func.SetParLimits(2, mean-100., mean+100.);
  func.SetParLimits(3, 2.,20.);
  Int_t ipar=ndp-1;
  if (_LET) {// LE tail
    func.SetParLimits(ipar+1, 0.05, 0.3);
    func.SetParLimits(ipar+2, 20., 50.);
    //func.FixParameter(ipar+1, fixed_tail_ratio);
    //func.FixParameter(ipar+2, fixed_tail_slope);
    ipar=ipar+2;
  }
  if (_HET) {// HE tail
    func.SetParLimits(ipar+1, 0.0, 0.2);
    func.SetParLimits(ipar+2, 2., 50.);
    ipar=ipar+2;
  }
  if (_BG) {// BG
    func.SetParLimits(ipar+1, 0., 1e5);
    func.FixParameter(ipar+2, 0.0);// for flat bg
    ipar=ipar+2;
  }

  // do fits
  h->Fit(funcname,"QLR0");// MIGRAD
  TString fitopt="ELR0";// Minos
  if (QUIET) fitopt="Q"+fitopt;// quiet mode
  h->Fit(funcname,fitopt);

  // fill function and params
  ff.SetFunction(func);
  Bool_t ENE=false;
  if (!QUIET) ff.PrintParams(ENE);
  if (PLOT) {
    TCanvas *cfit = ff.Plot(h,ENE);
    if (SAVE) {
      if (_fout) _fout->cd(_hname);
      cfit->Write();
      cfit->Close();
    }
  }
    
  h->Delete();
  return ff.GetFitParam();
}
  


void EneCalib::DoCalibE62(Bool_t PSSAVE=false,
                          Bool_t FITPLOT=false, Bool_t FITSAVE=false, Bool_t QUIET=false,
                          Bool_t CALSAVE=false){
  // peak positions or peaksearch
  if (_gmean) {
    Int_t gnp = _gmean->GetN();
    Double_t *p = _gmean->GetY();
    _pp.clear();
    for (int i=0; i<gnp; i++) {
      _pp.push_back(p[i]);
      //cout << Form("graph mean [%d]: %.1f",i,p[i]) << endl;
    }
  } else {
    DoPeakSearch(10,6,0.020,3,PSSAVE);
  }
  // calib lines and fit functions
  const Int_t np = int(_pp.size());
  Double_t ene[np],  eneErr[np];
  Double_t mean[np], meanErr[np];
  const Int_t np_zero = int(_pp.size())+1;
  Double_t ene_zero[np_zero],  eneErr_zero[np_zero];
  Double_t mean_zero[np_zero], meanErr_zero[np_zero];
  ene_zero[0]=0.;   eneErr_zero[0]=0.;
  mean_zero[0]=0.;  meanErr_zero[0]=0.;
  std::vector<fitparam> fps{};
  std::vector<TString> peaknames{};
  if (np==5 || np==6) {
    peaknames.push_back("CrKAlpha");
    peaknames.push_back("CrKBeta");
    peaknames.push_back("CoKAlpha");
    peaknames.push_back("CoKBeta");
    peaknames.push_back("CuKAlpha");
    if (np==6) peaknames.push_back("CuKBeta");
    _lines.clear();
    for (Int_t i=0; i<np; i++) {
      _lines.push_back(peaknames[i]);
    }
    SetFitFunctions(_lines);
  } else {
    cout << "number of peaks looks strange " << _hname << endl;
    ofstream ofst(badpsfile, ios::app);
    if (!ofst) {cout << "file open error " << badpsfile << endl; exit(1);}
    ofst << getdt() << " " << _hname << " " << np << endl;
    return;
  }
    
  for(Int_t i=0; i != np; ++i) {
    fps.push_back(DoFitRaw(_ffs[i],_pp[i],_lines[i].Data(),FITPLOT,FITSAVE,QUIET));
    mean[i]    = fps[i].mean;
    meanErr[i] = fps[i].meanErr;
    ene[i]     = fps[i].ene;
    eneErr[i]  = fps[i].eneErr;
    mean_zero[i+1]    = fps[i].mean;
    meanErr_zero[i+1] = fps[i].meanErr;
    ene_zero[i+1]     = fps[i].ene;
    eneErr_zero[i+1]  = fps[i].eneErr;
  }

  Clear();
  _gmean       = GetGraphErrors(np,ene,mean,eneErr,meanErr,_hname,gmean_tag);
  _sp3mean     = GetSpline3(_gmean,_hname,sp3mean_tag);
  _gsp3meander = GetSpline3DerivativeGraph(_sp3mean,_hname,1000,4000.,10000.);
  _gcalib      = GetGraphErrors(np,mean,ene,meanErr,eneErr,_hname,gcalib_tag);
  _sp3calib    = GetSpline3(_gcalib,_hname,sp3calib_tag);
  _hene        = GetEnergyConvertedHist(_hraw,_sp3calib,_hname,hene_tag);
    
  _gmean_zero       = GetGraphErrors(np_zero,ene_zero,mean_zero,eneErr_zero,meanErr_zero,_hname_zero,gmean_tag);
  _sp3mean_zero     = GetSpline3(_gmean_zero,_hname_zero,sp3mean_tag);
  _gsp3meander_zero = GetSpline3DerivativeGraph(_sp3mean_zero,_hname_zero,2000,0.,10000.);
  _gcalib_zero      = GetGraphErrors(np_zero,mean_zero,ene_zero,meanErr_zero,eneErr_zero,_hname_zero,gcalib_tag);
  _sp3calib_zero    = GetSpline3(_gcalib_zero,_hname_zero,sp3calib_tag);
  _hene_zero        = GetEnergyConvertedHist(_hraw,_sp3calib_zero,_hname_zero,hene_tag);

  if (CALSAVE) {
    if (_fout) _fout->cd(_hname);
    _hraw        ->Write();
    _gmean       ->Write();
    _sp3mean     ->Write();
    _gsp3meander ->Write();
    _gcalib      ->Write();
    _sp3calib    ->Write();
    _gmean_zero       ->Write();
    _sp3mean_zero     ->Write();
    _gsp3meander_zero ->Write();
    _gcalib_zero      ->Write();
    _sp3calib_zero    ->Write();
    _hene        ->Write();
    _hene_zero   ->Write();
    
    TString calibname=Form("%s_%s",calib_tag.Data(),_hname.Data());
    TCanvas *ccalib = new TCanvas(calibname,calibname,1,1,600,800);
    ccalib->Divide(1,2);
    ccalib->cd(1);
    gPad->SetGrid(0,0);
    gPad->SetTicks();
    SetAxisTGraphErrors(_gmean,"Energy (eV)","adc (ch)",61,61);
    _gmean->GetXaxis()->SetLimits(ene[0]-1000,ene[np-1]+1000);
    _gmean->SetMarkerStyle(20);
    _gmean->Draw("AP");
    _sp3mean->SetLineColor(4);
    _sp3mean->Draw("same");
    TLegend *lmean = new TLegend(0.15, 0.65, 0.5, 0.80);
    lmean->SetFillColor(0);
    lmean->SetFillStyle(0);
    lmean->SetTextSize(0.03);
    lmean->AddEntry(_gmean->GetName(),"fitted means","ep");
    lmean->AddEntry(_sp3mean->GetName(),"cubic spline interpolate","l");
    lmean->Draw("same");
    ccalib->cd(2);
    gPad->SetGrid(0,0);
    gPad->SetTicks();
    SetAxisTGraph(_gsp3meander,"Energy (eV)","eV/ch",61,61);
    _gsp3meander->GetXaxis()->SetLimits(ene[0]-1000,ene[np-1]+1000);
    _gsp3meander->SetMinimum(0.45);
    _gsp3meander->SetMaximum(0.90);
    _gsp3meander->SetLineColor(4);
    _gsp3meander->Draw("APC");
    TLegend *le2c = new TLegend(0.15, 0.65, 0.5, 0.80);
    le2c->SetFillColor(0);
    le2c->SetFillStyle(0);
    le2c->SetTextSize(0.04);
    le2c->AddEntry(_gsp3meander->GetName(),"e2c from 1./cubic_spline_derivative","l");
    le2c->Draw("same");
    ccalib->Update();
    //------------------------------------------
    TString calibname_zero=Form("%s_%s",calib_tag.Data(),_hname_zero.Data());
    TCanvas *ccalib_zero = new TCanvas(calibname_zero,calibname_zero,1,1,600,800);
    ccalib_zero->Divide(1,2);
    ccalib_zero->cd(1);
    gPad->SetGrid(0,0);
    gPad->SetTicks();
    SetAxisTGraphErrors(_gmean_zero,"Energy (eV)","adc (ch)",61,61);
    _gmean_zero->GetXaxis()->SetLimits(ene_zero[0]-1000,ene_zero[np_zero-1]+1000);
    _gmean_zero->SetMarkerStyle(20);
    _gmean_zero->Draw("AP");
    _sp3mean_zero->SetLineColor(4);
    _sp3mean_zero->Draw("same");
    TLegend *lmean_zero = new TLegend(0.15, 0.65, 0.5, 0.80);
    lmean_zero->SetFillColor(0);
    lmean_zero->SetFillStyle(0);
    lmean_zero->SetTextSize(0.03);
    lmean_zero->AddEntry(_gmean_zero->GetName(),"fitted means","ep");
    lmean_zero->AddEntry(_sp3mean_zero->GetName(),"cubic spline interpolate","l");
    lmean_zero->Draw("same");
    ccalib_zero->cd(2);
    gPad->SetGrid(0,0);
    gPad->SetTicks();
    SetAxisTGraph(_gsp3meander_zero,"Energy (eV)","eV/ch",61,61);
    _gsp3meander_zero->GetXaxis()->SetLimits(ene_zero[0]-1000,ene_zero[np_zero-1]+1000);
    _gsp3meander_zero->SetMinimum(0.45);
    _gsp3meander_zero->SetMaximum(0.90);
    _gsp3meander_zero->SetLineColor(4);
    _gsp3meander_zero->Draw("APC");
    TLegend *le2c_zero = new TLegend(0.15, 0.65, 0.5, 0.80);
    le2c_zero->SetFillColor(0);
    le2c_zero->SetFillStyle(0);
    le2c_zero->SetTextSize(0.04);
    le2c_zero->AddEntry(_gsp3meander_zero->GetName(),"e2c from 1./cubic_spline_derivative","l");
    le2c_zero->Draw("same");
    ccalib_zero->Update();
      
    ccalib->Write();
    ccalib_zero->Write();
    ccalib->Close();
    ccalib_zero->Close();
  }
    
  return;
}





#endif
