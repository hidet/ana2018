#ifndef CONVERT2HISTENE_CXX
#define CONVERT2HISTENE_CXX 1

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


inline TH1F* LoopChs(TString hntag="on",
                     TString hntag_calib="sprmon",
                     Int_t run=381,Bool_t ZERO=false,
                     TString fname="hoge.root",
                     TFile *fout=nullptr)
{
  if (gSystem->AccessPathName(fname)) {
    cout << "Error: cannot open file " << fname << endl;
    exit(1);
  }
  TString fcalib=Form("%s/%s/%s/%s_%s%d.root",
                      par_dir.Data(),calib_tag.Data(),hntag_calib.Data(),
                      calib_tag.Data(),hntag_calib.Data(),run);
  if (gSystem->AccessPathName(fcalib)) {
    cout << "Error: cannot open file " << fcalib << endl;
    exit(1);
  }
  // -----------------------------------------------------------
  TFile f(fname,"read");
  TList *list=f.GetListOfKeys();
  TFile fc(fcalib,"read");
  TList *listfc=fc.GetListOfKeys();
  f.cd();
  TObjLink *lnk=list->FirstLink();
  TString namecheck=Form("%s_%s%d_ch",hpht_phc.Data(),hntag.Data(),run);
  TList *hene_list=new TList();
  // -----------------------------------------------------------
  cout << "calib loop start" << endl;
  while (lnk) {
    TString keyname=lnk->GetObject()->GetName();
    if (keyname.Contains(namecheck)) {
      f.cd();
      TH1F *h1=(TH1F*)f.Get(keyname);
      TString hname=h1->GetName();
      Ssiz_t start=hname.Last('_');// <- the last '_' should be channel
      Ssiz_t end=hname.Length();
      TString chnum( hname(start+1,end) );
      TString hname_calib=Form("%s_%s%d_%s",hpht_phc.Data(),hntag_calib.Data(),run,chnum.Data());
      if (ZERO) hname_calib+="_zero";
      TString sp3name=Form("%s_%s",sp3calib_tag.Data(),hname_calib.Data());
      if (! listfc->Contains(sp3name) ) {
        cout << "cannot find " << sp3name << endl;
        if (h1) h1->Delete();
        lnk = lnk->Next();
        continue;
      }
      TSpline3 *sp3calib=(TSpline3*)fc.Get(sp3name);
      cout << Form("EneCalib: %s with %s",hname.Data(),hname_calib.Data()) << endl;
      TH1F *hene = EnergyConversion(h1,sp3calib);
      //TH1F *hene = EnergyConversion(h1,hname_calib,ZERO);
      if (!hene) cout << "something bad, not added " << chnum << endl;
      else if (hene->GetEntries()>0) hene_list->Add(hene);
      if (h1) h1->Delete();
      if (sp3calib) sp3calib->Delete();
    }
    lnk = lnk->Next();
  }
  fc.cd();
  fc.Close();
  cout << "calib loop end" << endl;
  // -----------------------------------------------------------
  TString hsum_name=Form("%s_%s%d_sum",hene_phc.Data(),hntag.Data(),run);
  if (ZERO) hsum_name += "_zero";
  TH1F* h0=(TH1F*)hene_list->First();
  TH1F* hene_sum=(TH1F*)h0->Clone(hsum_name);
  hene_sum->Reset();
  hene_sum->SetTitle(hsum_name);
  hene_sum->Merge(hene_list);
  cout << hene_sum->GetName() << " " << hene_sum->GetEntries() << endl;
  if (fout) {
    fout->cd();
    hene_sum->Write();
  }
  return hene_sum;
}


inline void LoopRuns(TString hntag_calib="sprmon",
                     TString target="hoge",
                     Bool_t ZERO=false)
{
  Int_t nruns=1;
  std::vector<int> runs;
  if (target=="he4") {
    nruns = nruns_he4;
    for (int i=0; i<nruns; ++i) runs.push_back(runs_he4[i]);
  } else if (target=="he3") {
    nruns = nruns_he3;
    for (int i=0; i<nruns; ++i) runs.push_back(runs_he3[i]);
  } else {// user define e.g.,
    runs={396,397};
    //runs={385,386,387,389};
    nruns=int(runs.size());
  }
  // ----------------------------------------------------
  TString foutname=Form("%s/%s_%s%d_%d_%s_spline",
                        dout.Data(),hene_phc.Data(),hntag_calib.Data(),
                        runs[0],runs[nruns-1],target.Data());
  if (ZERO) foutname+="_zero";
  foutname+=".root";
  cout << "runs: ";
  for (int i=0; i<nruns; ++i) cout << runs[i] << " "; cout << endl;
  cout << "foutname: " << foutname << endl;
  //check_continue();
  //backup_rootfile(foutname);
  TFile *fout = new TFile(foutname,"recreate");  
  // ----------------------------------------------------
  TString hntag="on";               TList *hene_on_sum_list = new TList();
  TString hntag_sprmon="sprmon";    TList *hene_sprmon_sum_list = new TList();
  TString hntag_sprmoff="sprmoff";  TList *hene_sprmoff_sum_list = new TList();
  TString hntag_kheton="kheton";    TList *hene_kheton_sum_list = new TList();
  TString hntag_khetoff="khetoff";  TList *hene_khetoff_sum_list = new TList();
  for (int i=0;i<nruns;++i) {
    TString frname=Form("%s/%s_%s%d_sum_spline",dout.Data(),hene_phc.Data(),hntag_calib.Data(),runs[i]);
    if (ZERO) frname+="_zero";
    frname+=".root";
    if (gSystem->AccessPathName(frname)) {
      cout << "missing file: " << frname << endl;
      continue;
    }
    TFile frun(frname,"read");
    fout->cd();
    cout << "reading: " << frname << endl;
    TString hn_on=Form("%s_%s%d_sum",hene_phc.Data(),hntag.Data(),runs[i]);
    if (ZERO) hn_on+="_zero";
    hene_on_sum_list->Add((TH1F*)frun.Get(hn_on)->Clone(Form("h_%s%d",hntag.Data(),runs[i])));
    TString hn_sprmon=Form("%s_%s%d_sum",hene_phc.Data(),hntag_sprmon.Data(),runs[i]);
    if (ZERO) hn_sprmon+="_zero";
    hene_sprmon_sum_list->Add((TH1F*)frun.Get(hn_sprmon)->Clone(Form("h_%s%d",hntag_sprmon.Data(),runs[i])));
    TString hn_sprmoff=Form("%s_%s%d_sum",hene_phc.Data(),hntag_sprmoff.Data(),runs[i]);
    if (ZERO) hn_sprmoff+="_zero";
    hene_sprmoff_sum_list->Add((TH1F*)frun.Get(hn_sprmoff)->Clone(Form("h_%s%d",hntag_sprmoff.Data(),runs[i])));
    TString hn_kheton=Form("%s_%s%d_sum",hene_phc.Data(),hntag_kheton.Data(),runs[i]);
    if (ZERO) hn_kheton+="_zero";
    hene_kheton_sum_list->Add((TH1F*)frun.Get(hn_kheton)->Clone(Form("h_%s%d",hntag_kheton.Data(),runs[i])));
    TString hn_khetoff=Form("%s_%s%d_sum",hene_phc.Data(),hntag_khetoff.Data(),runs[i]);
    if (ZERO) hn_khetoff+="_zero";
    hene_khetoff_sum_list->Add((TH1F*)frun.Get(hn_khetoff)->Clone(Form("h_%s%d",hntag_khetoff.Data(),runs[i])));
    frun.Close();
  }
  // ----------------------------------------------------
  fout->cd();
  TString hene_on_sum_name=Form("%s_%s%d_%d_runs",hene_phc.Data(),hntag.Data(),runs[0],runs[nruns-1]);
  if (ZERO) hene_on_sum_name += "_zero";
  TH1F* hon0=(TH1F*)hene_on_sum_list->First();
  TH1F* hene_on_sum=(TH1F*)hon0->Clone(hene_on_sum_name);
  hene_on_sum->Reset();
  hene_on_sum->SetTitle(hene_on_sum_name);
  hene_on_sum->Merge(hene_on_sum_list);
  fout->cd(); hene_on_sum->Write();
  // ----------------------------------------------------
  TString hene_sprmon_sum_name=Form("%s_%s%d_%d_runs",hene_phc.Data(),hntag_sprmon.Data(),runs[0],runs[nruns-1]);
  if (ZERO) hene_sprmon_sum_name += "_zero";
  TH1F* hon1=(TH1F*)hene_sprmon_sum_list->First();
  TH1F* hene_sprmon_sum=(TH1F*)hon1->Clone(hene_sprmon_sum_name);
  hene_sprmon_sum->Reset();
  hene_sprmon_sum->SetTitle(hene_sprmon_sum_name);
  hene_sprmon_sum->Merge(hene_sprmon_sum_list);
  fout->cd(); hene_sprmon_sum->Write();
  // ----------------------------------------------------
  TString hene_sprmoff_sum_name=Form("%s_%s%d_%d_runs",hene_phc.Data(),hntag_sprmoff.Data(),runs[0],runs[nruns-1]);
  if (ZERO) hene_sprmoff_sum_name += "_zero";
  TH1F* hon2=(TH1F*)hene_sprmoff_sum_list->First();
  TH1F* hene_sprmoff_sum=(TH1F*)hon2->Clone(hene_sprmoff_sum_name);
  hene_sprmoff_sum->Reset();
  hene_sprmoff_sum->SetTitle(hene_sprmoff_sum_name);
  hene_sprmoff_sum->Merge(hene_sprmoff_sum_list);
  fout->cd(); hene_sprmoff_sum->Write();
  // ----------------------------------------------------
  TString hene_kheton_sum_name=Form("%s_%s%d_%d_runs",hene_phc.Data(),hntag_kheton.Data(),runs[0],runs[nruns-1]);
  if (ZERO) hene_kheton_sum_name += "_zero";
  TH1F* hon3=(TH1F*)hene_kheton_sum_list->First();
  TH1F* hene_kheton_sum=(TH1F*)hon3->Clone(hene_kheton_sum_name);
  hene_kheton_sum->Reset();
  hene_kheton_sum->SetTitle(hene_kheton_sum_name);
  hene_kheton_sum->Merge(hene_kheton_sum_list);
  fout->cd(); hene_kheton_sum->Write();
  // ----------------------------------------------------
  TString hene_khetoff_sum_name=Form("%s_%s%d_%d_runs",hene_phc.Data(),hntag_khetoff.Data(),runs[0],runs[nruns-1]);
  if (ZERO) hene_khetoff_sum_name += "_zero";
  TH1F* hon4=(TH1F*)hene_khetoff_sum_list->First();
  TH1F* hene_khetoff_sum=(TH1F*)hon4->Clone(hene_khetoff_sum_name);
  hene_khetoff_sum->Reset();
  hene_khetoff_sum->SetTitle(hene_khetoff_sum_name);
  hene_khetoff_sum->Merge(hene_khetoff_sum_list);
  fout->cd(); hene_khetoff_sum->Write();
  // ----------------------------------------------------  
  cout << foutname << " has been created" << endl;  
  cout << "closed" << endl;
  fout->Close();
  return;
}




int convert2histene(TString hntag_calib="sprmon",
                    TString target="hoge",
                    Bool_t ZERO=false)
{
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
    runs={396,397};
    //runs={385,386,387,389};
    nruns=int(runs.size());
  }
  // ----------------------------------------------------
  if (gSystem->AccessPathName(fname)) {
    cout << "Error: cannot open file " << fname << endl;
    exit(1);
  }
  cout << "fname: " << fname << endl;
  cout << "runs: ";
  for (int i=0; i<nruns; ++i) cout << runs[i] << " "; cout << endl;
  check_continue();
  // ----------------------------------------------------
  TString hntag="on";             
  TString hntag_sprmon="sprmon";  
  TString hntag_sprmoff="sprmoff";
  TString hntag_kheton="kheton";  
  TString hntag_khetoff="khetoff";
  for (int i=0;i<nruns;++i) {
    TString foutname_1run=Form("%s/%s_%s%d_sum_spline",
                               dout.Data(),hene_phc.Data(),hntag_calib.Data(),runs[i]);
    if (ZERO) foutname_1run+="_zero";
    foutname_1run+=".root";
    backup_rootfile(foutname_1run);
    TFile fout1run(foutname_1run,"recreate");
    TH1F* hene_on      = LoopChs(hntag,hntag_calib,runs[i],ZERO,fname,&fout1run);
    TH1F *hene_sprmon  = LoopChs(hntag_sprmon,hntag_calib,runs[i],ZERO,fname,&fout1run);
    TH1F *hene_sprmoff = LoopChs(hntag_sprmoff,hntag_calib,runs[i],ZERO,fname,&fout1run);
    TH1F *hene_kheton  = LoopChs(hntag_kheton,hntag_calib,runs[i],ZERO,fname,&fout1run);
    TH1F *hene_khetoff = LoopChs(hntag_khetoff,hntag_calib,runs[i],ZERO,fname,&fout1run);
    fout1run.Close();
    cout << foutname_1run << " has been created" << endl;
  }
  LoopRuns(hntag_calib,target,ZERO);
  return 0;
}


#endif
