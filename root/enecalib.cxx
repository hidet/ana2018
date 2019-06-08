#ifndef ENERGYCALIB_CXX
#define ENERGYCALIB_CXX 1

#include <vector>
#include <iostream>

#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TDirectory.h"

#include "util.h"
#include "enecalib.hh"


inline void LoopChs(TString hntag="onsprmon", Int_t run=381,
                    TString fname="input.root", Int_t rebin=4,
                    Bool_t LET=true,Bool_t HET=false,Bool_t BG=false,
                    Bool_t SAVE=true, TString foutname="output.root")
{
  if (gSystem->AccessPathName(fname)) {
    cout << "Error: cannot open file " << fname << endl;
    exit(1);
  }
  TFile f(fname,"read");
  TDirectory *hist = (TDirectory*)f.GetDirectory("hist");
  TList *list = hist->GetListOfKeys();
  TObjLink *lnk = list->FirstLink();
  TString checkname=Form("%s_%s%d_ch",hpht_phc.Data(),hntag.Data(),run);
  Int_t ich=0;
  TString fprename=foutname;
  while (lnk) {
    TString keyname = lnk->GetObject()->GetName();
    if (keyname.Contains(checkname)) {
      if (ich==0) {
        fprename = backup_rootfile(foutname);// backup
        gSystem->Exec(Form("rm -f %s",foutname.Data()));
      }   
      TH1F *h1 = (TH1F*)hist->Get(keyname);
      EneCalib ec(h1,rebin,LET,HET,BG);
      if (fprename!=foutname) ec.SetPreviousCalib(fprename);// access the backup file
      if (SAVE) ec.OpenOutputFile(foutname);// update mode
      // PSSAVE, FITPLOT, FITSAVE, QUIET, CALSAVE
      ec.DoCalibE62(true,false,false,true,false);// peak search, rough calib
      ec.DoCalibE62(false,true,true,false,SAVE);
      if (SAVE) ec.CloseOutputFile();
      //TH1F *hene = ec.GetHistEne();
      //TH1F *hene_zero = ec.GetHistEneZero();
      h1->Delete();
      ich+=1;
    }
    lnk=lnk->Next();
  }
  f.cd();
  f.Close();
  return;
}


int enecalib(TString target="hoge", TString beam="on",
             TString spill="off", Int_t pre=0, Int_t post=500)
{
  gROOT->SetBatch(true);
  //---------------------------------------------------------
  if (beam!="on" && beam!="off") {
    cout << "beam should be off or on" << endl;
    return 1;
  }
  //---------------------------------------------------------
  Int_t   nruns=1;
  TString hntag="";
  TString spilltag="";
  std::vector<int> runs;
  Int_t rebin=4;
  Bool_t LET=true;Bool_t HET=false;Bool_t BG=false; Bool_t SAVE=true;
  //---------------------------------------------------------
  if (target=="he4") {// He4
    nruns=nruns_he4;
    for (int i=0; i<nruns; ++i) runs.push_back(runs_he4[i]);
  } else if (target=="he3") {// He3
    nruns=nruns_he3;
    for (int i=0; i<nruns; ++i) runs.push_back(runs_he3[i]);
  } else {// Test
    runs={160};
    nruns=int(runs.size());
  }
  if (beam=="on")        hntag="onsprmon";
  else if (beam=="off")  hntag="offsprmon";
  if (pre>0 || post>0) spilltag+=Form("_pre%03d_post%03d",pre,post);
  if (spill=="on")       spilltag+="_spillon";
  else if (spill=="off") spilltag+="_spilloff";
  spilltag+="_sprmcon_jbrscon_prime";// for E62 dump root
  
  cout << hntag << endl;
  cout << spilltag << endl;
  cout << "runs: ";
  for (int i=0; i<nruns; ++i) cout << runs[i] << " ";
  cout << endl;
  check_continue();

  TString fdir="", fname="", foutname="";
  for (int i=0;i<nruns;++i) {
    fdir  = Form("%s/run%04d",rootdir.Data(),runs[i]);
    void *dirp = gSystem->OpenDirectory(fdir);
    if (dirp) {
      TString fntmp="";
      while (fntmp!=nullptr) {
        fntmp = gSystem->GetDirEntry(dirp);
        if ( fntmp.Contains("mass") && fntmp.Contains(spilltag) ) {
          fname = Form("%s/%s",fdir.Data(),fntmp.Data());
          break;
        }
      }
      if (fntmp==nullptr) {
        cout << "Error: no mass file exists in " << fdir << endl;
        continue;
      }
    } else {
      cout << "Error: cannot open directory " << fdir << endl;
      continue;
    }
    foutname=Form("%s/run%04d_calib_%s%s.root",fdir.Data(),runs[i],hntag.Data(),spilltag.Data());
    cout << "input:  " << fname << endl;
    cout << "output: " << foutname << endl;
    LoopChs(hntag,runs[i],fname,rebin,LET,HET,BG,SAVE,foutname);
  }
  return 0;
}

#endif
