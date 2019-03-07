#ifndef FITUSER_HH
#define FITUSER_HH 1

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


class TF1;
class TH1;
struct xrayline;
struct fitparam;
class FitUser
{
public:
  FitUser(){}
  FitUser(const xrayline xl,
          const Bool_t let=false,const Bool_t het=false,const Bool_t bg=false,
          const TString funcname="",const TString tag="",
          const Double_t lr=0,const Double_t hr=30000,const Int_t npx=10000) {
    xl_=xl;
    npx_=npx;
    funcname_=funcname;
    tag_=tag;
    fp_.element=xl_.GetElement();
    fp_.linetype=xl_.GetLineType();    
    fp_.LET=let;
    fp_.HET=het;
    fp_.BG=bg;
    fp_.lrange=lr;
    fp_.hrange=hr;
    ndp_=5;// default number of parameters
    Int_t np=ndp_;
    Int_t np_bg=2;// <- need to change if you use other bg function
    if (fp_.LET) {np+=2;}
    if (fp_.HET) {np+=2;}
    if (fp_.BG)  {np+=np_bg;}
    fp_.npar=np;
  }
  ~FitUser(){}

  Double_t operator() (Double_t *x, Double_t *par)
  {
    // set parameters
    fp_.e2c=par[0];
    fp_.ampl=par[1];
    fp_.mean=par[2];
    fp_.gaussian_sigma=par[3];
    fp_.width=par[4];// <- this is new [eV]
    Int_t ipar=ndp_-1;// 5-1=4
    if (fp_.LET) {
      fp_.letail_ratio=par[ipar+1];
      fp_.letail_beta=par[ipar+2]/fp_.e2c;
      ipar=ipar+2;
    } else {
      fp_.letail_ratio=0;
      fp_.letail_beta=1;
    }
    if (fp_.HET) {
      fp_.hetail_ratio=par[ipar+1];
      fp_.hetail_beta=par[ipar+2]/fp_.e2c;
      ipar=ipar+2;
    } else {
      fp_.hetail_ratio=0;
      fp_.hetail_beta=1;
    }
    if (fp_.BG) {
      fp_.bg_0=par[ipar+1];
      fp_.bg_1=par[ipar+2];
      ipar=ipar+2;
    } else {
      fp_.bg_0=0;
      fp_.bg_1=0;
    }
    Double_t result=0;
    for(Int_t i=0; i != int(xl_.energies.size()); ++i) {
      Double_t mi = (xl_.energies[i]-xl_.energies[0])/fp_.e2c+fp_.mean;
      //Double_t wi = xl_.lorentzian_widths[i]/fp_.e2c;
      Double_t wi = fp_.width/fp_.e2c;// for free width
      Double_t ii = fp_.ampl*xl_.normalized_intensities[i]/xl_.normalized_intensities[0];
      if (fp_.LET && fp_.HET) {
        result += voigt(x,ii*(1.-fp_.letail_ratio-fp_.hetail_ratio),mi,fp_.gaussian_sigma,wi);
        result += letail(x,ii*fp_.letail_ratio,mi,fp_.gaussian_sigma,fp_.letail_beta);
        result += hetail(x,ii*fp_.hetail_ratio,mi,fp_.gaussian_sigma,fp_.hetail_beta);
      } else if (fp_.LET && !fp_.HET) {
        result += voigt(x,ii*(1.-fp_.letail_ratio),mi,fp_.gaussian_sigma,wi);
        result += letail(x,ii*fp_.letail_ratio,mi,fp_.gaussian_sigma,fp_.letail_beta);
      } else if (!fp_.LET && fp_.HET) {
        result += voigt(x,ii*(1.-fp_.hetail_ratio),mi,fp_.gaussian_sigma,wi);
        result += hetail(x,ii*fp_.hetail_ratio,mi,fp_.gaussian_sigma,fp_.hetail_beta);
      } else {
        result += voigt(x,ii,mi,fp_.gaussian_sigma,wi);
      }
    }
    if (fp_.BG) result += fp_.bg_0+fp_.bg_1*x[0];
    return result;
  }
  
  void SetFunction(TF1 &f1) {
    func_=TF1(f1);
    SetPrivateParams();
  }
  void SetFuncParams(Double_t *par, const Double_t *parErr) {
    func_.SetParameters(par);
    func_.SetParErrors(parErr);
  }
  void SetFuncParaNames(TString namePar[]) {
    for (Int_t i=0; i<fp_.npar; i++) {func_.SetParName(i,namePar[i].Data());}
  }
  void SetPrivateParams();
  void SetXrayLine(xrayline xl) {xl_=xl;}
  void SetFitParam(fitparam fp) {fp_=fp;}
  void SetNpar(Int_t np) {fp_.npar=np;}
  void SetFuncRange(Double_t lr, Double_t hr) {fp_.lrange=lr; fp_.hrange=hr;}
  void SetNpX(Int_t np) {npx_=np; func_.SetNpx(npx_);}
  void SetNpX() {func_.SetNpx(npx_);}
  
  TF1      GetFunction() {return func_;}  
  Int_t    GetNpar() {return fp_.npar;}  
  Double_t GetFuncLowRange() {return fp_.lrange;}
  Double_t GetFuncHighRange() {return fp_.hrange;}
  Int_t    GetNpX() {return npx_;}
  TString  GetParName(Int_t i) {return func_.GetParName(i);}
  xrayline GetXrayLine() {return xl_;}
  fitparam GetFitParam() {return fp_;}

  void Plot();
  void Plot(TH1F *h,Bool_t ENE,Bool_t SAVE);
  void Print();
  void PrintParams(Bool_t ENE){fp_.Print(ENE);}
  void SaveParams(Bool_t ENE);
  

private:
  xrayline xl_;// structure of xrayline (xrayline.hh)
  TString  tag_;// for save
  fitparam fp_;// structure of fit parameters (fitparam.hh)
  TF1      func_;// main fit function
  TString  funcname_;// name of function
  Int_t    npx_;// number of points of function
  Int_t    ndp_;// number of default parameters

};

void FitUser::SetPrivateParams()
{
  Double_t *par = func_.GetParameters();
  const Double_t *parErr = func_.GetParErrors();
  fp_.e2c   = par[0];
  fp_.e2cErr  = parErr[0];
  fp_.ampl  = par[1];
  fp_.amplErr = parErr[1];
  fp_.mean  = par[2];
  fp_.meanErr = parErr[2];
  fp_.gaussian_sigma=par[3];
  fp_.gaussian_sigmaErr = parErr[3];
  fp_.width = par[4];
  fp_.widthErr = parErr[4];
  fp_.fwhm = fp_.gaussian_sigma*fp_.e2c*2.355;
  fp_.fwhmErr = GetError(fp_.gaussian_sigma,fp_.gaussian_sigmaErr,
                         fp_.e2c,fp_.e2cErr,"MUL")*2.355;
  fp_.ene = xl_.energies[0];
  fp_.eneErr = xl_.energies_err[0];
  fp_.enediff = fp_.mean-xl_.energies[0];
  fp_.enediffErr = GetError(fp_.mean,fp_.meanErr,
                            xl_.energies[0],xl_.energies_err[0],"SUB");
  fp_.chisq = func_.GetChisquare();
  fp_.ndf = (Double_t)func_.GetNDF();
  fp_.nchi = fp_.chisq/fp_.ndf;
  Int_t ipar=ndp_-1;// 5-1=4
  if (fp_.LET) {
    fp_.letail_ratio  = par[ipar+1];
    fp_.letail_ratioErr = parErr[ipar+1];
    fp_.letail_beta   = par[ipar+2];
    fp_.letail_betaErr  = parErr[ipar+2];
    ipar=ipar+2;
  }
  if (fp_.HET) {
    fp_.hetail_ratio=par[ipar+1];
    fp_.hetail_ratioErr = parErr[ipar+1];
    fp_.hetail_beta   = par[ipar+2];
    fp_.hetail_betaErr  = parErr[ipar+2];
    ipar=ipar+2;
  }
  if (fp_.BG) {
    fp_.bg_0  = par[ipar+1];
    fp_.bg_0Err = parErr[ipar+1];
    fp_.bg_1  = par[ipar+2];
    fp_.bg_1Err = parErr[ipar+2];
    ipar=ipar+2;
  }
}

void FitUser::Print(){
  Double_t *par = func_.GetParameters();
  const Double_t *parErr = func_.GetParErrors();
  cout << endl;
  cout << Form("### Fit ###") << endl;
  cout << Form("  npar = %d",fp_.npar) << endl;
  for (Int_t i=0; i<fp_.npar; i++) {
    cout << Form("  par[%d] %s = %f +- %f",i,func_.GetParName(i),par[i],parErr[i]) << endl;
  }
}

void FitUser::Plot(TH1F *h,Bool_t ENE=false, Bool_t SAVE=true)
{
  Double_t bw=h->GetBinWidth(1);
  TString caname = Form("c_%s_%s",h->GetName(),tag_.Data());
  TCanvas *cfit = new TCanvas(caname,caname,1,1,801,601);
  gPad->SetGrid(0,0);
  gPad->SetTicks();
  if (!ENE) {SetAxisHist(h,"adc [ch]",Form("Counts / %.2f ch",bw),61,61);}
  else {SetAxisHist(h,"Energy [eV]",Form("Counts / %.2f eV",bw),61,61);}
  Double_t maxh = h->GetMaximum();
  h->GetYaxis()->SetRangeUser(0.,maxh*1.5);
  h->DrawCopy("e");
  func_.SetLineWidth(2);
  func_.SetLineColor(3);
  func_.DrawCopy("same");
  TLegend *lhist = new TLegend(0.58, 0.45, 0.85, 0.80);
  lhist->SetFillColor(0);
  lhist->SetFillStyle(0);
  lhist->SetTextSize(0.03);
  lhist->AddEntry(h->GetName(),h->GetName(),"pe");
  lhist->AddEntry(func_.GetName(),Form("%s",tag_.Data()),"l");
  lhist->AddEntry((TObject*)0,Form("FWHM  = %.2f+-%.2f (eV)",fp_.fwhm,fp_.fwhmErr),"");
  lhist->AddEntry((TObject*)0,Form("width = %.2f+-%.2f (eV)",fp_.width,fp_.widthErr),"");
  if (!ENE) {
    lhist->AddEntry((TObject*)0,Form(" e2c = %.2f+-%.2f (eV/ch)",fp_.e2c,fp_.e2cErr),"");
    lhist->AddEntry((TObject*)0,Form("mean = %.2f+-%.2f (ch)",fp_.mean,fp_.meanErr),"");
  }else {
    lhist->AddEntry((TObject*)0,Form("mean = %.2f+-%.2f (eV)",fp_.mean,fp_.meanErr),"");
    lhist->AddEntry((TObject*)0,Form("diff = %.4f+-%.4f (eV)",fp_.enediff,fp_.enediffErr),"");
  }
  if (fp_.LET) {
    lhist->AddEntry((TObject*)0,Form("letail %% %.2f+-%.2f ",fp_.letail_ratio*100,fp_.letail_ratioErr*100),"");
    lhist->AddEntry((TObject*)0,Form("letail beta %.2f+-%.2f ",fp_.letail_beta,fp_.letail_betaErr),"");
  }
  if (fp_.HET) {
    lhist->AddEntry((TObject*)0,Form("hetail %% %.2f+-%.2f ",fp_.hetail_ratio*100,fp_.hetail_ratioErr*100),"");
    lhist->AddEntry((TObject*)0,Form("hetail beta %.2f+-%.2f ",fp_.hetail_beta,fp_.hetail_betaErr),"");
  }
  if (fp_.BG) {
    lhist->AddEntry((TObject*)0,Form("bg const %.2f+-%.2f ",fp_.bg_0,fp_.bg_0Err),"");
    lhist->AddEntry((TObject*)0,Form("bg slope %.2f+-%.2f ",fp_.bg_1,fp_.bg_1Err),"");
  }
  //lhist->SetNColumns(2);
  lhist->Draw("same");
  cfit->Update();
  if (SAVE) {
    TString coutname=Form("%s/%s/%s_%s",fig_dir.Data(),fitUser_tag.Data(),funcname_.Data(),tag_.Data());
    if (ENE) {coutname+="_ene";}
    coutname+=".pdf";
    cfit->SaveAs(coutname);
    //lhist->Delete();
    //cfit->Close();
  }
}

void FitUser::SaveParams(Bool_t ENE=false)
{
  TString foutname=Form("%s/%s/%s_%s",par_dir.Data(),fitUser_tag.Data(),funcname_.Data(),tag_.Data());
  if (ENE) {foutname+="_ene";}
  foutname+=".dat";
    ofstream file(foutname, ios::out | ios::binary);
  if (!file) {
    cerr << "Error: cannot open file " << foutname << endl;
    exit(1);
  }else {
    file.write((char*)&fp_, sizeof(fp_));
    cout << foutname << " has been created" << endl;
  }
}



inline fitparam FitUserXray(TH1F *h1,
                            Int_t nrebin=1,Double_t pp=0,TString ln="KHe4LAlpha",
                            Bool_t LET=true,Bool_t HET=false,Bool_t BG=false,
                            Bool_t ENE=false,Bool_t PLOT=true,Bool_t SAVE=true)
                           
{
  // h1: pointer of histogram, refer should work
  // nrebin: rebin to fit for the histogram
  // pp: peak position, if ENE, pp is not used
  // ln: name of the xray line (e.g., MnKAlpha, MnKa, mnka)
  // LET: bool for low-energy tail on/off
  // HET: bool for high-energy tail on/off
  // BG: bool for background function on/off
  // ENE: bool, true if the x-axis of h1 is energy

  // use a copy to keep the original histogram
  TString hname = h1->GetName();// <- this is VERY important
  TString htitle = h1->GetTitle();
  TH1F *h = (TH1F*)h1->Clone(Form("%s_f",hname.Data()));
  h->Rebin(nrebin);
  h->SetTitle(htitle);
  Double_t bw=h->GetBinWidth(1);
  
  // xrayline
  xrayline xl = get_xl(ln);
  //xl.Print();
  TString element = xl.GetElement();
  TString linetype = xl.GetLineType();

  // parameters
  Double_t ampl=1000.;
  Double_t e2c=1.0;
  Double_t mean=xl.energies[0];
  Double_t width=xl.lorentzian_widths[0];
  if (!ENE) {
    e2c=0.6;
    mean=pp;// from peakserch()
  }
  Double_t gaussian_sigma = 5;
    
  // initial parameter arrays
  std::vector<string> paranames{"eV per ch","Voigt ampl","Voigt mean","Voigt sigma","Voigt gamma"};
  std::vector<double> tp{e2c,ampl,mean,gaussian_sigma,width};
  Int_t ndp=Int_t(tp.size());
  std::vector<double> tpErr(ndp,0.);// 5 elements with 0.
  if (LET) {
    paranames.push_back("LE tail ratio");
    paranames.push_back("LE tail beta");
    tp.push_back(0.1);// LE tail ratio
    tp.push_back(50.);// LE tail beta
    tpErr.push_back(0.);// LE tail ratio
    tpErr.push_back(0.);// LE tail beta
  }
  if (HET) {
    paranames.push_back("HE tail ratio");
    paranames.push_back("HE tail beta");
    tp.push_back(0.1);// HE tail ratio
    tp.push_back(50.);// HE tail beta
    tpErr.push_back(0.);// HE tail ratio
    tpErr.push_back(0.);// HE tail beta
  }
  if (BG) {
    paranames.push_back("BG const");
    paranames.push_back("BG slope");
    tp.push_back(0.);// BG const
    tp.push_back(0.);// BG slope
    tpErr.push_back(0.);// BG const
    tpErr.push_back(0.);// BG slope
  }
  // range
  h->SetAxisRange(mean*0.987,mean*1.013,"X");
  Double_t low = mean*0.99;
  Double_t high = mean*1.01;
  // fixed parameters
  Double_t fixed_sigma=6.5/2.355;
  Double_t fixed_tail_ratio=0.20;
  Double_t fixed_tail_slope=30.;
  if ( (ln=="KHe4LAlpha") || (ln=="KHe4La")) {
    low=mean*0.997;
    high=mean*1.005;
    fixed_sigma=6.5/2.355;
    fixed_tail_ratio=0.19;
    fixed_tail_slope=20.;
  } else if ( (ln=="KHe3LAlpha")  || (ln=="KHe3La")) {
    low=mean*0.997;
    high=mean*1.005;
    fixed_sigma=7.2/2.355;
    fixed_tail_ratio=0.19;
    fixed_tail_slope=27.;
  }
  // fitfunc
  Int_t npx=10000;// number of drawing points of function
  TString funcname = Form("%s_%s",fitUser_tag.Data(),hname.Data());
  TString xl_tag=Form("%s%s",element.Data(),linetype.Data());
  FitUser ff(xl,LET,HET,BG,funcname,xl_tag,low,high,npx);
  Int_t npar = ff.GetNpar();
  TF1 func(funcname,ff,low,high,npar);
  func.SetNpx(npx);
  func.SetParameters(&tp.front());
  func.SetParErrors(&tpErr.front());
  for (Int_t i=0; i<npar; i++) {func.SetParName(i,paranames[i].data());}
  // fit parameter boundary
  if (!ENE) {
    func.SetParLimits(0, e2c*0.5, e2c*1.5);
    if (linetype=="KBeta") func.FixParameter(0, e2c);
  } else {func.FixParameter(0, 1.0);}
  func.SetParLimits(1, 10., 1e8);
  func.SetParLimits(2, mean-100., mean+100.);
  //func.SetParLimits(3, 2.,20.);
  func.FixParameter(3, fixed_sigma);// sigma
  func.SetParLimits(4, 0.,10.);// width
  Int_t ipar=ndp-1;// 5-1=4
  if (LET) {// LE tail
    //func.SetParLimits(ipar+1, 0.05, 0.3);
    //func.SetParLimits(ipar+2, 20., 50.);
    func.FixParameter(ipar+1, fixed_tail_ratio);
    func.FixParameter(ipar+2, fixed_tail_slope);
    ipar=ipar+2;
  }
  if (HET) {// HE tail
    func.SetParLimits(ipar+1, 0.0, 0.2);
    func.SetParLimits(ipar+2, 5., 50.);
    ipar=ipar+2;
  }
  if (BG) {// BG
    func.SetParLimits(ipar+1, -1e5, 1e5);
    //func.SetParLimits(ipar+2, -1e3, 1e3);
    func.FixParameter(ipar+2, 0.0);// for flat bg
    ipar=ipar+2;
  }

  // do fits
  h->Fit(funcname,"QLR0");
  h->Fit(funcname,"ELR0");

  // print, plot, and save fit results
  ff.SetFunction(func);
  ff.PrintParams(ENE);
  if (SAVE) {ff.SaveParams(ENE);}
  if (PLOT) {ff.Plot(h,ENE,SAVE);}

  h->Delete();
  
  return ff.GetFitParam();
}



#endif
