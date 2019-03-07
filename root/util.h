#ifndef UTIL_H
#define UTIL_H 1

#include <iostream>
#include <string>
#include <time.h>

#include "TString.h"
#include "TSystem.h"

using namespace std;

// --------------------------------------------------------
// He3
const Int_t nruns_he3=96;
Int_t runs_he3[nruns_he3] = {160,161,162,163,165,166,167,170,171,172,
                             175,176,177,178,179,180,181,182,183,186,
                             187,188,189,190,191,192,193,194,195,196,
                             197,198,201,202,203,204,205,206,207,208,
                             209,210,212,215,216,217,218,222,223,224,
                             225,226,227,228,231,232,233,239,240,241,
                             242,259,260,261,263,264,266,267,268,269,
                             270,271,272,275,276,277,278,279,281,282,
                             283,284,286,287,290,291,292,293,294,295,
                             296,297,298,299,300,301};
// --------------------------------------------------------
// He4
const Int_t nruns_he4=70;
Int_t runs_he4[nruns_he4] = {320,321,327,328,329,330,331,332,333,334,
                             335,336,337,338,340,341,345,346,347,348,
                             349,350,351,352,353,354,355,356,999,359,
                             363,364,366,367,368,369,370,371,372,373,
                             374,375,376,377,381,382,383,384,385,386,
                             387,389,395,396,397,398,399,400,401,402,
                             403,406,407,408,410,411,421,422,423,424};
// --------------------------------------------------------

// save directories
static const TString par_dir="./par";
static const TString fig_dir="./fig";

static const TString badpsfile="./memo/bad_peack_search.txt";

std::vector<string> badch{"ch257","ch265","ch293"};

static const TString dout="./output";
static const TString makechhist_dir="pyroot";
static const TString anadir=gSystem->Getenv("HEATESANADIR");
static const TString fhe3_on=Form("%s/%s/output/hpht_phc_on_run0160_0301.root",anadir.Data(),makechhist_dir.Data());
static const TString fhe4_on=Form("%s/%s/output/hpht_phc_on_run0320_0424.root",anadir.Data(),makechhist_dir.Data());
static const TString fhe3_off=Form("%s/%s/output/hpht_phc_off_run0160_0301.root",anadir.Data(),makechhist_dir.Data());
static const TString fhe4_off=Form("%s/%s/output/hpht_phc_off_run0320_0424.root",anadir.Data(),makechhist_dir.Data());

// important tags used for directories and file names
static const TString hpht_phc       = "hpht_phc";
static const TString hene_phc       = "hene_phc";
static const TString fitTES_tag     = "fitTES";// fit function
static const TString fitElem_tag    = "fitTES";// fit function
static const TString fitUser_tag    = "fitTES";// fit function
static const TString ps_tag         = "peak_search";// peak search
static const TString sp3meander_tag = "sp3meander";// energy calibration
static const TString sp3mean_tag    = "sp3mean";// energy calibration
static const TString gmean_tag      = "gmean";// energy calibration
static const TString gcalib_tag     = "gcalib";// energy calibration
static const TString sp3calib_tag   = "sp3calib";// energy calibration
static const TString hene_tag       = "hene";// energy converted histo
static const TString calib_tag      = "calib";// output root file


// C++ datetime and timestamp
// from http://yut.hatenablog.com/entry/20130617/1371425713
#define NT 30
inline string getdt() {
  // get datetime now
   time_t timer = time(NULL);
   struct tm* tm = localtime(&timer);
   char datetime[NT];
   strftime(datetime,NT,"%Y-%m-%d %H:%M:%S",tm);
   return datetime;
}

inline string getdt2() {
  // get datetime now
   time_t timer = time(NULL);
   struct tm* tm = localtime(&timer);
   char datetime[NT];
   strftime(datetime,NT,"%Y%m%d%H%M%S",tm);
   return datetime;
}

inline unsigned int getts() {
  // get timestamp now
  time_t epoch_time = time(NULL);
  return epoch_time;
}

inline unsigned int dt2ts(const char* datetime="2013-06-15 12:00:00") {
  // datetime to timestamp
  struct tm tm;
  if( strptime(datetime,"%Y-%m-%d %H:%M:%S",&tm) != NULL ) {
    return mktime(&tm);
  }
  return 0;
}

inline string ts2dt(unsigned int timestamp=1371265200) {
  // timestamp to datetime
  time_t timer=(time_t)timestamp;
  struct tm* tm=localtime(&timer);
  char datetime[NT];
  strftime(datetime,NT,"%Y-%m-%d %H:%M:%S",tm);
  return datetime;
}

inline Bool_t contains_check(TString test,std::vector<string> names)
{
  Bool_t flag=false;
  for (int i=0; i<(int)names.size(); i++) {
    if ( test.Contains(names[i]) ) flag=true;
  }
  return flag;
}

inline Bool_t backup_rootfile(TString fname)
{
  if (!gSystem->AccessPathName(fname)) {
    Ssiz_t end = fname.Last('.');// get the position of extention
    TString fnamebkp( fname(0,end) );
    fnamebkp+=Form("_%s.root",getdt2().data());
    gSystem->Exec(Form("cp -p %s %s",fname.Data(),fnamebkp.Data()));
    cout << Form("backup file: %s has been created",fnamebkp.Data()) << endl;
    return true;
  }
  return false;
}

inline Bool_t check_continue()
{
  Int_t loop=0;
  while(loop==0){
    cout << "This will take a long time. Do you want to do? (y/n)" << endl;
    char a; cin >> a;
    if (a=='y' || a=='Y'){
      loop=1;
      return true;
    }
    else if(a=='n' || a=='N'){
      loop=1;
      exit(1);
      return false;
    }
    else{
      cout <<"please answer y or n" << endl;
      loop=0;
    }
  }
  return false;
}


#endif
